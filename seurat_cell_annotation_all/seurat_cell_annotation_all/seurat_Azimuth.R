#https://app.azimuth.hubmapconsortium.org/app/human-pbmc 

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

# * Cell annotations (at multiple levels of resolution)
# * Prediction scores (i.e. confidence scores) for each annotation
# * Projection onto the reference-derived 2-dimensional UMAP visualization

#OLD OBJECT (with subset function)
mtx_cell <- readRDS("/suffolk/WorkGenomicsD/ec474/mitox_cellcycle/mtx_cellcycle_regressed_clustered_atac_all_singlets_mgatk.rds")

#Seurat object
DefaultAssay(mtx_cell) <- "RNA"


#Rename WNN outputs
mtx_cell@meta.data$predicted.celltype.l1.wnn <- mtx_cell@meta.data$predicted.celltype.l1
mtx_cell@meta.data$predicted.celltype.l2.wnn <- mtx_cell@meta.data$predicted.celltype.l2
mtx_cell@meta.data$predicted.celltype.l1.score.wnn <- mtx_cell@meta.data$predicted.celltype.l1.score
mtx_cell@meta.data$predicted.celltype.l2.score.wnn <- mtx_cell@meta.data$predicted.celltype.l2.score


# The RunAzimuth function can take a Seurat object as input
mtx_cell <- RunAzimuth(mtx_cell, reference = "pbmcref")

# Warning: Overwriting miscellanous data for model
# Warning: Adding a dimensional reduction (refUMAP) without the associated assay being present
# Warning: Adding a dimensional reduction (refUMAP) without the associated assay being present
# detected inputs from HUMAN with id type Gene.name
# reference rownames detected HUMAN with id type Gene.name
# Normalizing query using reference SCT model
# Warning: 113 features of the features specified were not present in both the reference query assays.
# Continuing with remaining 4887 features.
# Projecting cell embeddings
# Finding query neighbors
# Finding neighborhoods
# Finding anchors
# 	Found 531 anchors
# Finding integration vectors
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Predicting cell labels
# Predicting cell labels
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Predicting cell labels
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#   |                                                  | 0 % ~calculating
# Integrating dataset 2 with reference dataset
# Finding integration vectors
# Integrating data
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s
# Computing nearest neighbors
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# Running UMAP projection
# 16:40:41 Read 401 rows
# 16:40:41 Processing block 1 of 1
# 16:40:41 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 20
# 16:40:41 Initializing by weighted average of neighbor coordinates using 1 thread
# 16:40:41 Commencing optimization for 67 epochs, with 8020 positive edges
# Using method 'umap'
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 16:40:41 Finished
# Projecting reference PCA onto query
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Projecting back the query cells into original PCA space
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Computing scores:
#     Finding neighbors of original query cells
#     Finding neighbors of transformed query cells
#     Computing query SNN
#     Determining bandwidth and computing transition probabilities
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Total elapsed time: 0.598118305206299

# Warning messages:
# 1: Different cells and/or features from existing assay prediction.score.celltype.l2
# 2: Keys should be one or more alphanumeric characters followed by an underscore, setting key from integrated_dr_ to integrateddr_
# 3: In RunUMAP.default(object = neighborlist, reduction.model = reduction.model,  :
#   Number of neighbors between query and reference is not equal to the number of neighbors within reference
# 4: No assay specified, setting assay as RNA by default.

#Rename Azimuth output
mtx_cell@meta.data$predicted.celltype.l1.azimuth <- mtx_cell@meta.data$predicted.celltype.l1
mtx_cell@meta.data$predicted.celltype.l2.azimuth <- mtx_cell@meta.data$predicted.celltype.l2
mtx_cell@meta.data$predicted.celltype.l3.azimuth <- mtx_cell@meta.data$predicted.celltype.l3
mtx_cell@meta.data$predicted.celltype.l1.score.azimuth <- mtx_cell@meta.data$predicted.celltype.l1.score
mtx_cell@meta.data$predicted.celltype.l2.score.azimuth <- mtx_cell@meta.data$predicted.celltype.l2.score
mtx_cell@meta.data$predicted.celltype.l3.score.azimuth <- mtx_cell@meta.data$predicted.celltype.l3.score

colnames(mtx_cell@meta.data)
# [71] "wsnn_res.0.8"                        "predicted.celltype.l1.score"
# [73] "predicted.celltype.l1"               "predicted.celltype.l2.score"
# [75] "predicted.celltype.l2"               "predicted.celltype.l1.wnn"
# [77] "predicted.celltype.l2.wnn"           "predicted.celltype.l1.score.wnn"
# [79] "predicted.celltype.l2.score.wnn"     "predicted.celltype.l3.score"
# [81] "predicted.celltype.l3"               "mapping.score"
# [83] "predicted.celltype.l1.azimuth"       "predicted.celltype.l2.azimuth"
# [85] "predicted.celltype.l3.azimuth"       "predicted.celltype.l1.score.azimuth"
# [87] "predicted.celltype.l2.score.azimuth" "predicted.celltype.l3.score.azimuth"

mtx_cell
# An object of class Seurat
# 145382 features across 401 samples within 7 assays
# Active assay: RNA (36601 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  6 other assays present: ATAC, HTO, SCT, prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
#  10 dimensional reductions calculated: pca, tsne, pca.cellcycle, pca.aftercell, umap.aftercell, umap.aftercell_15, lsi, umap.atac, integrated_dr, ref.umap

p1 <- DimPlot(mtx_cell, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3)
ggsave(p1, filename="reference_celltypel1_azimuth.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_azimuth/")

p2 <- DimPlot(mtx_cell, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
ggsave(p2, filename="reference_celltypel2_azimuth.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_azimuth/")

p3 <- DimPlot(mtx_cell, group.by = "predicted.celltype.l3", label = TRUE, label.size = 3) 
ggsave(p2, filename="reference_celltypel3_azimuth.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_azimuth/")


vln <- VlnPlot(mtx_cell, group.by = "predicted.celltype.l1", features="predicted.celltype.l1.score")
ggsave(vln, filename="vln_celltypel1_azimuth.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_azimuth/")


vln <- VlnPlot(mtx_cell, group.by = "predicted.celltype.l2", features="predicted.celltype.l2.score")
ggsave(vln, filename="vln_celltypel2_azimuth.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_cell_annotation_all/figures_azimuth/")


table(mtx_cell@meta.data$predicted.celltype.l1.azimuth, mtx_cell@meta.data$predicted.celltype.l1.wnn)
#             B CD4 T CD8 T  DC Mono  NK other T
#   B        71     0     0   0    0   0       0
#   CD4 T     3   278     4   0    0   0       0
#   CD8 T     0    19    79   0    0  10       1
#   DC        0     0     0   3    0   0       0
#   Mono      0     0     0   0    2   0       0
#   NK        0     0     1   0    0 123       0
#   other T   0     7     2   0    0   0      31

#                ASDC B intermediate B memory B naive CD16 Mono CD4 Naive
#   B intermediate    0              9        7      40         0         0
#   B naive           0              0        0      14         0         0
#   CD16 Mono         0              0        0       0         2         0
#   CD4 Naive         0              0        0       0         0        12
#   CD4 TCM           0              0        0       4         0        56
#   CD8 Naive         0              0        0       0         0         0
#   CD8 TCM           0              0        0       0         0         0
#   CD8 TEM           0              0        0       1         0         1
#   cDC2              1              0        0       0         0         0
#   MAIT              0              0        0       0         0         0
#   NK                0              0        0       0         0         0
#   NK_CD56bright     0              0        0       0         0         0
#   pDC               0              0        0       0         0         0

#                  CD4 TCM CD4 TEM CD8 Naive CD8 TCM CD8 TEM dnT gdT MAIT  NK
#   B intermediate       0       0         0       0       0   0   0    0   0
#   B naive              0       0         0       0       0   0   0    0   0
#   CD16 Mono            0       0         0       0       0   0   0    0   0
#   CD4 Naive            0       0         1       0       0   0   0    0   0
#   CD4 TCM            189      20         2       0       2   1   0    0   0
#   CD8 Naive            0       0        11       0       0   0   0    0   0
#   CD8 TCM              1       1         2       4       2   0   0    0   0
#   CD8 TEM              2       6         0       1      67   0   1    0  10
#   cDC2                 0       0         0       0       0   0   0    0   0
#   MAIT                 3       4         0       0       3   0   0   31   0
#   NK                   0       0         0       0       0   0   0    0 112
#   NK_CD56bright        0       0         0       0       1   0   0    0   7
#   pDC                  0       0         0       0       0   0   0    0   0

#                  NK_CD56bright pDC
#   B intermediate             0   0
#   B naive                    0   0
#   CD16 Mono                  0   0
#   CD4 Naive                  0   0
#   CD4 TCM                    0   0
#   CD8 Naive                  0   0
#   CD8 TCM                    0   0
#   CD8 TEM                    0   0
#   cDC2                       0   0
#   MAIT                       0   0
#   NK                         0   0
#   NK_CD56bright              1   0
#   pDC                        0   2

table_l2 <- table(mtx_cell@meta.data$predicted.celltype.l2.azimuth, mtx_cell@meta.data$predicted.celltype.l2.wnn)
write.csv(table_l2, "/suffolk/WorkGenomicsD/ec474/mitox_azimuth/table_l2_azimuth_wnn.csv")

# We can visualize the expression of canonical marker genes to examine the accuracy of predictions. Note that Azimuth normalizes data (internally) before mapping, but does not return the results, so we normalize the data here before visualization.

# Here, we specifically visualize:

# * The expression of CCR7 on CD4 and CD8 Naive T cells
# * The expression of FCGR3A on CD16+ monocytes, CD56dim NK cells, and cytotoxic CD8 T cells
# * The expression of AXL on rare populations of AXL+SIGLEC6+ dendritic cells ([ASDC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5775029/))
# * Prediction scores for the annotation CD4+ regulatory T cells (Treg)

mtx_cell <- NormalizeData(mtx_cell)
Idents(mtx_cell) <- "predicted.celltype.l2"

p1 <- FeaturePlot(mtx_cell, features = "CCR7")

p2 <- FeaturePlot(mtx_cell, features = "FCGR3A")

p3 <- VlnPlot(mtx_cell, features = "AXL", group.by = "predicted.celltype.l2", idents = c("ASDC", "pDC", "cDC1", "cDC2"))
#None of the cells requested found in this assay

p4 <- FeaturePlot(mtx_cell, features = "predictionscorecelltypel2_Treg")

p1 + p2 + p4 + plot_layout(ncol=3)
ggsave(filename="features_azimuth.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_WNN/figures_azimuth/")

saveRDS(mtx_cell, file="/suffolk/WorkGenomicsD/ec474/mitox_azimuth/mtx_singlets_mgatk_cellcycleregressed_annotated_wnn_azimuth.rds")


