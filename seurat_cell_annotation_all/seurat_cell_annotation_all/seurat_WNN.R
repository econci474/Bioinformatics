#load libraries
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(clustree)
library(patchwork)


#load object 
mtx_cell <- readRDS("/suffolk/WorkGenomicsD/ec474/mitox_cellcycle/mtx_cellcycle_regressed_clustered_atac_all_singlets_mgatk.rds")

DefaultAssay(mtx_cell) <- "RNA"

mtx_cell
# An object of class Seurat
# 279346 features across 634 samples within 5 assays
# Active assay: RNA (36601 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  4 other assays present: ATAC, HTO, mito, SCT
#  10 dimensional reductions calculated: pca.genes, pca.cellcycle, pca.aftercell, umap.aftercell, 
#  umap.aftercell_05, umap.aftercell_07, umap.aftercell_15, umap.aftercell_2, lsi, umap.atac

length(rownames(mtx_cell)) #genes
#[1] 36601

head(rownames(mtx_cell))
# "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" "AL627309.2"

gene_names <- rownames(mtx_cell)
gene_names[grep("_", gene_names)]
#character(0)

length(colnames(mtx_cell)) #cells
#[1] 254
# "AAACGTACACAGGGAC-1" "AAAGCACCACTCAACA-1" "AAAGGCTCAATCTCTC-1"

#Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. 
#We use this graph for UMAP visualization and clustering

mtx_cell <- FindMultiModalNeighbors(mtx_cell, reduction.list = list("pca.aftercell", "lsi"), dims.list = list(1:50, 2:50))
mtx_cell <- RunUMAP(mtx_cell, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mtx_cell <- FindClusters(mtx_cell, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

plot <- DimPlot(mtx_cell, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN") 
ggsave(plot, filename="umap_wnn.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_wnn/")

#Annotate the clusters below.
#Note that you could also annotate the dataset using our supervised mapping pipelines, 
#using either our [vignette](multimodal_reference_mapping.html), https://satijalab.org/seurat/articles/covid_sctmapping or 
#[automated web tool, Azimuth](https://satijalab.org/azimuth).

#Multimodal reference mapping using PBMC Cite-seq multimodal reference: https://zenodo.org/records/7779017#.ZCMojezMJqs 
reference <- readRDS("/suffolk/RawGenomicsF/ec474/multiome_reference/pbmc_multimodal_2023.rds")

# An object of class Seurat
# 20957 features across 161764 samples within 2 assays
# Active assay: SCT (20729 features, 3000 variable features)
#  2 layers present: counts, data
#  1 other assay present: ADT
#  6 dimensional reductions calculated: apca, aumap, pca, spca, umap, wnn.umap

head(rownames(reference)) #gene symbols
# "AL627309.1" "AL669831.5" "LINC00115"  "FAM41C"     "NOC2L" "KLHL17"

head(colnames(reference)) #What does the L1 in front of it mean? Cells 
# [1] "L1_AAACCCAAGAAACTCA" "L1_AAACCCAAGACATACA" "L1_AAACCCACAACTGGTT"
# [4] "L1_AAACCCACACGTACTA" "L1_AAACCCACAGCATACT" "L1_AAACCCACATCAGTCA"
# [99988] "E2L3_TCGGGCACAGGACGAT" "E2L3_TCGGGCAGTCTATGAC" "E2L3_TCGGGCAGTGATACTC"
# [99991] "E2L3_TCGGGCAGTGCCGAAA" "E2L3_TCGGGCAGTGTCACAT" "E2L3_TCGGGCAGTTCCTTGC"
# [99994] "E2L3_TCGGGCAGTTGGACTT" "E2L3_TCGGGCATCATTCACT" "E2L3_TCGGGCATCCCATACC"
# [99997] "E2L3_TCGGGCATCTCTAAGG" "E2L3_TCGGGCATCTTCGGTC" "E2L3_TCGGGTGAGGTCCCTG"

DefaultAssay(mtx_cell) <- "SCT"

mtx_cell
# An object of class Seurat
# 279346 features across 634 samples within 5 assays
# Active assay: SCT (11983 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  4 other assays present: RNA, ATAC, HTO, mito
#  11 dimensional reductions calculated: pca.genes, pca.cellcycle, pca.aftercell, umap.aftercell, umap.aftercell_05, 
#  umap.aftercell_07, umap.aftercell_15, umap.aftercell_2, lsi, umap.atac, wnn.umap

anchor <- FindTransferAnchors(
  reference = reference,
  query = mtx_cell,
  reference.reduction = "spca",
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT",
  dims = 1:50
)

#Normalizing query using reference SCT model
#Projecting cell embeddings
#Finding neighborhoods
#Finding anchors
#	Found 635 anchors

DefaultAssay(mtx_cell) <- "SCT"

mtx_cell <- MapQuery(
  anchorset = anchor,
  query = mtx_cell,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"
  ),
  reduction.model = "wnn.umap"
)

# Finding integration vectors
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Predicting cell labels
# Warning: Layer counts isn't present in the assay object; returning NULL
# Predicting cell labels
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Layer counts isn't present in the assay object; returning NULL
# Warning: Layer counts isn't present in the assay object; returning NULL
# Warning: Layer counts isn't present in the assay object; returning NULL
#   |                                                  | 0 % ~calculating
# Integrating dataset 2 with reference dataset
# Finding integration vectors
# Integrating data
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=07s
# Computing nearest neighbors
# Running UMAP projection
# 16:06:05 Read 634 rows
# 16:06:05 Processing block 1 of 1
# 16:06:05 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 15
# 16:06:05 Initializing by weighted average of neighbor coordinates using 1 thread
# 16:06:05 Commencing optimization for 67 epochs, with 9510 positive edges
# Using method 'umap'
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 16:06:05 Finished




#Cell annotation plots 
plot <- DimPlot(mtx_cell, reduction = "wnn.umap", group.by = "predicted.celltype.l1", alpha = 1, label = TRUE, label.size = 3)
ggsave(plot, filename="reference_celltypel1_wnn.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_wnn/")


plot <- DimPlot(mtx_cell, reduction = "wnn.umap", group.by = "predicted.celltype.l2", alpha = 1, label = TRUE, label.size = 3)
ggsave(plot, filename="reference_celltypel2_wnn.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_wnn/")


vln <- VlnPlot(mtx_cell, group.by = "predicted.celltype.l1", features="predicted.celltype.l1.score")
ggsave(vln, filename="vln_celltypel1_wnn.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_wnn/")


vln <- VlnPlot(mtx_cell, group.by = "predicted.celltype.l2", features="predicted.celltype.l2.score")
ggsave(vln, filename="vln_celltypel2_wnn.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_wnn/")





#Visualise and contrast umap.aftercell_15, umap.atac, wnn.umap
p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


#Do cell markers FeaturePlot (lavel by continuous variable) and Vln plots 
#DimPlot (lael by categorical feature)

#DEG of cell vs another cell within control group
