#Repeating workflow using Seurat object that has had RNA and ATAC QC done to it beforehand. 

library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)

#HTO expression matrix
data_dir <- "/suffolk/RawGenomicsF/ec474/cite-seq_count4/umi_count"
hto <- Read10X(data.dir = data_dir, gene.column = 1, cell.column = 1)
#gene.colum is 1 as feature.tsv.gz has format, and we want the hash ID, not the sequence

#zcat modified_features.tsv.gz | head
#NPC_HASH1	TAGTAATGGCAGGAA
#NPC_HASH2	AACAATGAGGTACGG
#NPC_HASH3	CTGAGGACATCTACC
#unmapped	NA

#cell.colum is 1 as barcode.tsv.gz only has one column 
#zcat barcodes.tsv.gz | head -n 10
#GGTTGGTGTGACATAT
#GGTCCGTAGCCTTAAA
#TCAAGTATCCCGAAGC
#GCGCTAGGTTTAGCGA
#GGCTGTCAGCCTGTGA
#ATGTTTGAGGGACCTC
#TATGACTCAGGCTTGT
#TCCCTCACACCTCACC
#TGCTTGCTCGATTATG
#TCACTGACACCTGGTG

#zcat matrix.mtx.gz | head -n 10
#%%MatrixMarket matrix coordinate integer general
#%
#4 10081 40317
#1 1 89
#3 1 61
#2 1 345
#4 1 1288
#2 2 248
#1 2 67
#3 2 112

dim(hto)
#[1]     4 10081

rowSums(hto)
#NPC_HASH1 NPC_HASH2 NPC_HASH3  unmapped
#   726,899   3,121,435   1,235,953  12,465,344

head(colnames(hto))
#[1] "GGTTATATCCCTCGCA" "CGGGACAAGCACAGAA" "CATAATGTCTAAATCG" "TGGGCCTAGCATCCAG"
#[5] "CTGGATGTCCAGGGAG" "TTTGTTGGTGATCAGC"

colnames(hto) <- paste0(colnames(hto), "-1")
# [1] "GGTTATATCCCTCGCA-1" "CGGGACAAGCACAGAA-1" "CATAATGTCTAAATCG-1"
# [4] "TGGGCCTAGCATCCAG-1" "CTGGATGTCCAGGGAG-1" "TTTGTTGGTGATCAGC-1"

mtx_umis <- readRDS("/suffolk/WorkGenomicsD/ec474/mitox_qc/mtx_rna_atac_qc.rds")

#An object of class Seurat
#134807 features across 9907 samples within 2 assays
#Active assay: RNA (36601 features, 0 variable features)
# 1 layer present: counts
# 1 other assay present: ATAC

head(colnames(mtx_umis))
#[1] "AAACAGCCACCTAATG-1" "AAACAGCCATAGACCC-1" "AAACAGCCATTAAAGG-1" "AAACATGCAAGGTGGC-1"
#[5] "AAACATGCAATCCTGA-1" "AAACATGCACCTCACC-1"

head(rownames(mtx_umis))
#[1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"
#[6] "AL627309.2"

#select barcodes present in both mtx RNA assay and hashtag matrix
joint_bcs <- intersect(colnames(mtx_umis), colnames(hto))
length(joint_bcs)
#[1] 9262 length of mtx_umis post quality controls

mtx_umis <- mtx_umis[, joint_bcs] 
hto <- as.matrix(hto[, joint_bcs])

# Normalize RNA data with log normalization
DefaultAssay(object = mtx_umis) <- "RNA"
mtx_umis <- NormalizeData(mtx_umis)
mtx_umis <- FindVariableFeatures(mtx_umis, selection.method = "mean.var.plot")
mtx_umis <- ScaleData(mtx_umis, features = VariableFeatures(mtx_umis))

#An object of class Seurat
#134807 features across 9907 samples within 2 assays
#Active assay: RNA (36601 features, 679 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: ATAC

#Add HTO matrix to multimodal Seurat object
mtx_umis[["HTO"]] <- CreateAssayObject(counts = hto)
#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

rownames(hto)
#[1] "NPC_HASH1" "NPC_HASH2" "NPC_HASH3" "unmapped"

hto_counts <- GetAssayData(mtx_umis, assay = "HTO", layer = "counts")
total_counts_per_cell <- colSums(hto_counts)

#head(total_counts_per_cell)
#AAACAGCCACCTAATG AAACAGCCATAGACCC AAACAGCCATTAAAGG AAACATGCAAGGTGGC
#            1204              968             1039             1187
#AAACATGCAATCCTGA AAACATGCACCTCACC
#            1101             1280

cells_to_keep <- names(total_counts_per_cell)[total_counts_per_cell > 0]
length(cells_to_keep)
#[1] 9262

#HTODemux with CLR normalization 
DefaultAssay(object = mtx_umis) <- "HTO"
mtx_umis <- NormalizeData(mtx_umis, assay = 'HTO', normalization.method = 'CLR')

str(mtx_umis)
#$ NormalizeData.HTO       :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#  .. .. .. ..@ name       : chr "NormalizeData.HTO"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2024-05-15 10:45:40"
#  .. .. .. ..@ assay.used : chr "HTO"
#  .. .. .. ..@ call.string: chr "NormalizeData(mtx_umis, assay = \"HTO\", normalization.method = \"CLR\")"
#  .. .. .. ..@ params     :List of 5
#  .. .. .. .. ..$ assay               : chr "HTO"
#  .. .. .. .. ..$ normalization.method: chr "CLR"
#  .. .. .. .. ..$ scale.factor        : num 10000
#  .. .. .. .. ..$ margin              : num 1
#  .. .. .. .. ..$ verbose             : logi TRUE

#https://satijalab.org/seurat/reference/htodemux 
#HTODemux(
#  object,
#  assay = "HTO",
#  positive.quantile = 0.99,
#  init = NULL,
#  nstarts = 100,
#  kfunc = "clara",
#  nsamples = 100,
#  seed = 42,
#  verbose = TRUE
#)

#HTODemux with k-means
mtx_umis_kmeans <- HTODemux(mtx_umis, assay = "HTO", positive.quantile = 0.99, init=4, kfunc="kmeans", nstarts=100)

#As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis: 
#https://satijalab.org/seurat/reference/aggregateexpression

#This message is displayed once per session.
#First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#This message is displayed once every 8 hours.

#Cutoff for NPC-HASH1 : 98 reads
#Cutoff for NPC-HASH2 : 431 reads
#Cutoff for NPC-HASH3 : 176 reads
#Cutoff for unmapped : 1071 reads

table(mtx_umis_kmeans$HTO_classification.global)
#Doublet Negative  Singlet
#    1209     5273     2780

#HTODemus with clara and k-medoids
mtx_umis_clara0.99_100 <- HTODemux(mtx_umis, assay = "HTO", positive.quantile = 0.99, init = 4, kfunc="clara", 
nsamples=100)
#First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#This message is displayed once every 8 hours.
#Cutoff for NPC-HASH1 : 98 reads
#Cutoff for NPC-HASH2 : 429 reads
#Cutoff for NPC-HASH3 : 175 reads
#Cutoff for unmapped : 1030 reads

table(mtx_umis_clara0.99_100$HTO_classification.global)
# Doublet Negative  Singlet
#    1214     5136     2912
#Total barcodes = 9907, which is the same as number of cells we expect

#I will stick with clara method since it gives me the greatest number of singlets. 
#Now I will try to optimise the number of singlets that can be obtained and cells tagegd by HASH-3 for my experiment

mtx_umis_clara_0.99_1000 <- HTODemux(mtx_umis, assay = "HTO", positive.quantile = 0.99, init = 4, kfunc="clara", 
nsamples=1000)
#Cutoff for NPC-HASH1 : 98 reads
#Cutoff for NPC-HASH2 : 431 reads
#Cutoff for NPC-HASH3 : 176 reads
#Cutoff for unmapped : 977 reads

table(mtx_umis_clara_0.99_1000$HTO_classification.global)
#Doublet Negative  Singlet
#    1212     4967     3083

table(mtx_umis_clara_0.99_1000$HTO_maxID)
#NPC-HASH1 NPC-HASH2 NPC-HASH3  unmapped
#     1746      1656      2377      3483

Idents(mtx_umis_clara_0.99_1000) <- 'HTO_classification.global'
VlnPlot(mtx_umis_clara_0.99_1000, features = 'nCount_RNA', pt.size = 0.01, log = TRUE)
ggsave(filename="vlnplot_clara0.99_1000_singlet_doublet_nCountRNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")

Idents(mtx_umis_clara_0.99_1000) <- 'HTO_maxID'
VlnPlot(mtx_umis_clara_0.99_1000, features = 'nCount_RNA', pt.size = 0.01, log = TRUE)
ggsave(filename="vlnplot_clara0.99_1000_hash_nCountRNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")


mtx_umis_clara_0.95_1000 <- HTODemux(mtx_umis, assay = "HTO", positive.quantile = 0.95, init = 4, kfunc="clara", 
nsamples=1000)
#Cutoff for NPC-HASH1 : 87 reads
#Cutoff for NPC-HASH2 : 375 reads
#Cutoff for NPC-HASH3 : 153 reads
#Cutoff for unmapped : 844 reads

table(mtx_umis_clara_0.95_1000$HTO_classification.global)
#Doublet Negative  Singlet
#    1212     4967     3083

table(mtx_umis_clara_0.95_1000$HTO_maxID)
#NPC-HASH1 NPC-HASH2 NPC-HASH3  unmapped
#     1746      1656      2377      3483

mtx_umis_clara_0.95_1000$HTO_classification.global <- factor(mtx_umis_clara_0.95_1000$HTO_classification.global, levels = c("Negative", "Singlet", "Doublet"))
Idents(mtx_umis_clara_0.95_1000) <- 'HTO_classification.global'
VlnPlot(mtx_umis_clara_0.95_1000, features = 'nCount_RNA', pt.size = 0.01, log = TRUE)
ggsave(filename="vlnplot_clara0.95_1000_singlet_doublet_nCountRNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")

Idents(mtx_umis_clara_0.95_1000) <- 'HTO_maxID'
VlnPlot(mtx_umis_clara_0.95_1000, features = 'nCount_RNA', pt.size = 0.01, log = TRUE)
ggsave(filename="vlnplot_clara0.95_1000_hash_nCountRNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")


mtx_umis_clara_0.90_1000 <- HTODemux(mtx_umis, assay = "HTO", positive.quantile = 0.90, init = 4, kfunc="clara", 
nsamples=1000)
#Cutoff for NPC-HASH1 : 81 reads
#Cutoff for NPC-HASH2 : 348 reads
#Cutoff for NPC-HASH3 : 142 reads
#Cutoff for unmapped : 779 reads

table(mtx_umis_clara_0.90_1000$HTO_classification.global)
#Doublet Negative  Singlet
#   2504     3586     3172

table(mtx_umis_clara_0.90_1000$HTO_maxID)
#NPC-HASH1 NPC-HASH2 NPC-HASH3  unmapped
#     1746      1656      2377      3483

Idents(mtx_umis_clara_0.90_1000) <- 'HTO_classification.global'
VlnPlot(mtx_umis_clara_0.90_1000, features = 'nCount_RNA', pt.size = 0.01, log = TRUE)
ggsave(filename="vlnplot_clara0.90_1000_singlet_doublet_nCountRNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")

Idents(mtx_umis_clara_0.90_1000) <- 'HTO_maxID'
VlnPlot(mtx_umis_clara_0.90_1000, features = 'nCount_RNA', pt.size = 0.01, log = TRUE)
ggsave(filename="vlnplot_clara0.90_1000_hash_nCountRNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")


#Plot 0.95 per 1000 group
# Group cells based on the max HTO signal
Idents(mtx_umis_clara_0.95_1000) <- 'HTO_maxID'
RidgePlot(mtx_umis_clara_0.95_1000, assay = 'HTO', features = rownames(mtx_umis[['HTO']])[1:4], ncol = 2)
ggsave(filename="ridgeplot_clara0.90_1000.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")

#Visualise hto scatter to ensure mutual exclusivity of singlets 
FeatureScatter(mtx_umis_clara_0.95_1000, feature1 = 'NPC-HASH1', feature2 = 'NPC-HASH2')
ggsave(filename="featurescatter_clara0.95_1000_hash1_hash2.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")

FeatureScatter(mtx_umis_clara_0.95_1000, feature1 = 'NPC-HASH2', feature2 = 'NPC-HASH3')
ggsave(filename="featurescatter_clara0.95_1000_hash2_hash3.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")


#Pick favourite HTODemux output and save it in RDS 
saveRDS(mtx_umis_clara_0.95_1000, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000.rds")


#Save the others too for later comparison
saveRDS(mtx_umis_clara_0.90_1000, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.90_1000.rds")
saveRDS(mtx_umis_clara_0.99_1000 , file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.99_1000.rds")



mtx_umis_clara_0.95_1000 <- readRDS("/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000.rds")
# An object of class Seurat
# 134811 features across 9262 samples within 3 assays
# Active assay: HTO (4 features, 0 variable features)
#  2 layers present: counts, data
#  2 other assays present: RNA, ATAC

#Subset away the negative cells
mtx_umis_clara_0.95_1000$HTO_classification.global <- factor(mtx_umis_clara_0.95_1000$HTO_classification.global, levels = c("Negative", "Doublet", "Singlet"))
DefaultAssay(mtx_umis_clara_0.95_1000) <- "HTO"
Idents(mtx_umis_clara_0.95_1000) <- 'HTO_classification.global'

unique(Idents(mtx_umis_clara_0.95_1000))
#Negative Doublet  Singlet
#Levels: Negative Doublet Singlet

table(mtx_umis_clara_0.95_1000$HTO_classification.global)
#Negative  Doublet  Singlet
#    4275     1799     3188

mtx_umis_clara_0.95_1000_subset <- subset(mtx_umis_clara_0.95_1000, idents = "Negative", invert = TRUE)

DefaultAssay(mtx_umis_clara_0.95_1000_subset) <- "HTO"
mtx_umis_clara_0.95_1000_subset <- ScaleData(mtx_umis_clara_0.95_1000_subset, features = rownames(mtx_umis_clara_0.95_1000_subset), verbose = FALSE)
mtx_umis_clara_0.95_1000_subset <- RunPCA(mtx_umis_clara_0.95_1000_subset, features = rownames(mtx_umis_clara_0.95_1000_subset), approx = FALSE)
mtx_umis_clara_0.95_1000_subset <- RunTSNE(mtx_umis_clara_0.95_1000_subset, dims = 1:4, perplexity = 100)

Idents(mtx_umis_clara_0.95_1000_subset) <- 'HTO_classification.global' #singlet vs douplet
DimPlot(mtx_umis_clara_0.90_1000_subset)
ggsave(filename="pca_singlet_douplet_clara0.95_1000.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")

Idents(mtx_umis_clara_0.95_1000_subset) <- 'HTO_maxID' #hashes
DimPlot(mtx_umis_clara_0.90_1000_subset)
ggsave(filename="pca_hash_clara0.95_1000.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")



#Subset the singlets
#Cluster and visualize cells using the usual scRNA-seq workflow, 
#and examine for the potential presence of batch effects.

# Extract the singlets
Idents(mtx_umis_clara_0.95_1000) <- 'HTO_classification.global'
mtx_umis_clara_0.95_1000_singlets <- subset(mtx_umis_clara_0.95_1000, idents = 'Singlet')
saveRDS(mtx_umis_clara_0.95_1000_singlets, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlets.rds")

mtx_umis_clara_0.95_1000_singlets <- readRDS(file="/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlets.rds")
Idents(mtx_umis_clara_0.95_1000_singlets) <- 'HTO_maxID'
mtx_umis_clara_0.95_1000_singlets_no_unmapped <- subset(mtx_umis_clara_0.95_1000_singlets, idents = 'unmapped', invert=TRUE)
saveRDS(mtx_umis_clara_0.95_1000_singlets_no_unmapped, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlets_no_unmapped.rds")


#mtx_umis_clara_0.90_1000_singlets
#An object of class Seurat
#134811 features across 3988 samples within 3 assays
#Active assay: RNA (36601 features, 679 variable features)
# 3 layers present: counts, data, scale.data
# 2 other assays present: ATAC, HTO

Idents(mtx_umis_clara_0.95_1000_singlets) <- 'HTO_maxID'
pdf("/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures/vlnplot_clara0.90_1000_singlets_hashtags_nCountRNA.pdf")
VlnPlot(mtx_umis_clara_0.95_1000_singlets, features = 'nCount_RNA', pt.size = 0.1, log = TRUE)
ggsave(filename="pca_clara0.95_1000_hash.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")


# Select the top 1000 most variable features
DefaultAssay(mtx_umis_clara_0.95_1000_singlets) <- "RNA"
mtx_umis_clara_0.95_1000_singlets <- FindVariableFeatures(mtx_umis_clara_0.95_1000_singlets, selection.method = 'mean.var.plot')

# Scaling RNA data, we only scale the variable features here for efficiency
mtx_umis_clara_0.95_1000_singlets <- ScaleData(mtx_umis_clara_0.95_1000_singlets, features = VariableFeatures(mtx_umis_clara_0.95_1000_singlets))

# Run PCA
mtx_umis_clara_0.95_1000_singlets <- RunPCA(mtx_umis_clara_0.95_1000_singlets, features = VariableFeatures(mtx_umis_clara_0.95_1000_singlets), reduction.name="pca_hto", reduction.key="htoPCA_")

ElbowPlot(mtx_umis_clara_0.95_1000_singlets, reduction="pca_hto")

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
mtx_umis_clara_0.95_1000_singlets <- FindNeighbors(mtx_umis_clara_0.95_1000_singlets, reduction = 'pca', dims = 1:10)
mtx_umis_clara_0.95_1000_singlets <- FindClusters(mtx_umis_clara_0.95_1000_singlets, resolution = 0.6, verbose = FALSE)
mtx_umis_clara_0.95_1000_singlets <- RunTSNE(mtx_umis_clara_0.95_1000_singlets, reduction = 'pca', dims = 1:10)

# Projecting singlet identities on TSNE visualization
pdf("/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures/PCA_clara0.90_1000_singlets.pdf")
DimPlot(mtx_umis_clara_0.90_1000_singlets, group.by = "HTO_classification")
ggsave(filename="pca_clara0.95_1000_singlets_doublets.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_hto_rna_atac")


#Save only NPC_HASH3 cells (irrespective of Negative/Singlet/Doublet)
Idents(mtx_umis_clara_0.95_1000) <- 'HTO_maxID'
mtx_umis_clara_0.95_1000_HASH3 <- subset(mtx_umis_clara_0.95_1000, idents = "NPC-HASH3")
saveRDS(mtx_umis_clara_0.95_1000_HASH3, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_HASH3.rds")

mtx_umis_clara_0.95_1000_HASH3


#Save only the NPC-HASH3 singlet cells
Idents(mtx_umis_clara_0.95_1000_singlets) <- 'HTO_maxID'
mtx_umis_clara_0.95_1000_singlets_HASH3 <- subset(mtx_umis_clara_0.95_1000_singlets, idents = "NPC-HASH3")
saveRDS(mtx_umis_clara_0.95_1000_singlets_HASH3, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlet_HASH3_RNA_ATAC_QC.rds")

mtx_umis_clara_0.95_1000_singlets_HASH3
#An object of class Seurat
#134811 features across 401 samples within 3 assays
#Active assay: RNA (36601 features, 1500 variable features)
# 3 layers present: counts, data, scale.data
# 2 other assays present: ATAC, HTO
# 2 dimensional reductions calculated: pca, tsne


#Save only the NPC-HASH1 singlet cells
Idents(mtx_umis_clara_0.95_1000_singlets) <- 'HTO_maxID'
mtx_umis_clara_0.95_1000_singlets_HASH1 <- subset(mtx_umis_clara_0.95_1000_singlets, idents = "NPC-HASH1")
saveRDS(mtx_umis_clara_0.95_1000_singlets_HASH1, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlet_HASH1_RNA_ATAC_QC.rds")

#Save only the NPC-HASH2 singlet cells
Idents(mtx_umis_clara_0.95_1000_singlets) <- 'HTO_maxID'
mtx_umis_clara_0.95_1000_singlets_HASH2 <- subset(mtx_umis_clara_0.95_1000_singlets, idents = "NPC-HASH2")
saveRDS(mtx_umis_clara_0.95_1000_singlets_HASH2, file = "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlet_HASH2_RNA_ATAC_QC.rds")
