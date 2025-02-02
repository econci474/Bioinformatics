library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(clustree)
library(patchwork)

#Lead in seurat object
mtx_cell <- readRDS(file="/suffolk/WorkGenomicsD/ec474/mitox_cellcycle/mtx_cellcycle_regressed_clustered_all_singlets_mgatk.rds")

#An object of class Seurat
#144167 features across 254 samples within 3 assays
#Active assay: SCT (9360 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 2 other assays present: RNA, ATAC
# 4 dimensional reductions calculated: pca.genes, pca.aftercell, umap.aftercell, umap.aftercell_15

#Check all unnecessary reductions deleted
mtx_cell@reductions

#metadata 
colnames(mtx_cell[[]])
# [1] "orig.ident"                        "nCount_RNA"
# [3] "nFeature_RNA"                      "nCount_ATAC"
# [5] "nFeature_ATAC"                     "percent.mt"
# [7] "nucleosome_signal"                 "nucleosome_percentile"
# [9] "TSS.enrichment"                    "TSS.percentile"
#[11] "gex_barcode"                       "atac_barcode"
#[13] "is_cell"                           "excluded_reason"
#[15] "gex_raw_reads"                     "gex_mapped_reads"
#[17] "gex_conf_intergenic_reads"         "gex_conf_exonic_reads"
#[19] "gex_conf_intronic_reads"           "gex_conf_exonic_unique_reads"
#[21] "gex_conf_exonic_antisense_reads"   "gex_conf_exonic_dup_reads"
#[23] "gex_exonic_umis"                   "gex_conf_intronic_unique_reads"
#[25] "gex_conf_intronic_antisense_reads" "gex_conf_intronic_dup_reads"
#[27] "gex_intronic_umis"                 "gex_conf_txomic_unique_reads"
#[29] "gex_umis_count"                    "gex_genes_count"
#[31] "atac_raw_reads"                    "atac_unmapped_reads"
#[33] "atac_lowmapq"                      "atac_dup_reads"
#[35] "atac_chimeric_reads"               "atac_mitochondrial_reads"
#[37] "atac_fragments"                    "atac_TSS_fragments"
#[39] "atac_peak_region_fragments"        "atac_peak_region_cutsites"
#[41] "pct_reads_in_peaks"                "blacklist_fraction"
#[43] "high.tss"                          "nucleosome_group"
#[45] "nCount_HTO"                        "nFeature_HTO"
#[47] "HTO_maxID"                         "HTO_secondID"
#[49] "HTO_margin"                        "HTO_classification"
#[51] "HTO_classification.global"         "hash.ID"
#[53] "nCount_SCT"                        "nFeature_SCT"
#[55] "S.Score"                           "G2M.Score"
# [57] "Phase"                             "old.ident"
# [59] "SCT_snn_res.1"                     "seurat_clusters"
# [61] "SCT_snn_res.0.2"                   "SCT_snn_res.0.5"
# [63] "SCT_snn_res.0.7"                   "SCT_snn_res.1.5"
# [65] "SCT_snn_res.2"

head(colnames(mtx_cell))
# [1] "AAACGCGCAGCACGTT-1" "AAACGCGCAGGCAAGC-1" "AAACGTACACAGGGAC-1"
# [4] "AAAGCACCACTCAACA-1" "AAAGGCTCAATCTCTC-1" "AAAGGTTAGTCTGGGC-1"

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mtx_cell) <- "ATAC"
mtx_cell <- FindTopFeatures(mtx_cell, min.cutoff = 'q0')
mtx_cell <- RunTFIDF(mtx_cell) #term frequency inverse document frequency (TF-IDF) normalization on a matrix
mtx_cell <- RunSVD(mtx_cell) #partial singular value decomposition using irlba
mtx_cell <- RunUMAP(mtx_cell, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Running SVD
# Scaling cell embeddings
# > mtx_cell <- RunUMAP(mtx_cell, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# 14:51:21 UMAP embedding parameters a = 0.9922 b = 1.112
# 14:51:21 Read 254 rows and found 49 numeric columns
# 14:51:21 Using Annoy for neighbor search, n_neighbors = 30
# 14:51:21 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 14:51:21 Writing NN index file to temp file /tmp/RtmpLkCBVS/file6cfd30a1c9a3
# 14:51:21 Searching Annoy index using 1 thread, search_k = 3000
# 14:51:21 Annoy recall = 100%
# 14:51:23 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 14:51:26 Initializing from normalized Laplacian + noise (using RSpectra)
# 14:51:26 Commencing optimization for 500 epochs, with 10186 positive edges
# Using method 'umap'
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 14:51:28 Optimization finished

plot <- DimPlot(mtx_cell, dims=c(2,3), reduction="lsi", label=TRUE) + labs(title="UMAP of ATAC")
ggsave(plot, filename="umap_atac_pca2v3.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_atac/")

plot <- DimPlot(mtx_cell, dims=c(3,4), reduction="lsi", label=TRUE) + labs(title="UMAP of ATAC")
ggsave(plot, filename="umap_atac_pca3v4.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_atac/")

plot <- DimPlot(mtx_cell, dims=c(4,5), reduction="lsi", label=TRUE) + labs(title="UMAP of ATAC")
ggsave(plot, filename="umap_atac_pca4v5.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_atac/")

#15 clusters found on ATAC, best umap is dims=c(3,4)

#An object of class Seurat
#144167 features across 254 samples within 3 assays
#Active assay: ATAC (98206 features, 98206 variable features)
# 2 layers present: counts, data
# 2 other assays present: RNA, SCT
# 6 dimensional reductions calculated: pca.genes, pca.aftercell, umap.aftercell, umap.aftercell_15, lsi, umap.atac

saveRDS(mtx_cell, file="/suffolk/WorkGenomicsD/ec474/mitox_cellcycle/mtx_cellcycle_regressed_clustered_atac_all_singlets_mgatk.rds")


#TF enrichment analysis - ATAC assay
#Scenic - RNA dat
#Scenic Plut - multiome RNA and ATAC

