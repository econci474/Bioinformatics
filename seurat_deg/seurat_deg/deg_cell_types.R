# The bulk of Seurat’s differential expression features can be accessed through the 
#FindMarkers() function. 

#By default, Seurat performs differential expression (DE) testing based on the 
#non-parametric Wilcoxon rank sum test. 

library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(patchwork)

#To test for DE genes between two specific groups of cells, 
#specify the ident.1 and ident.2 parameters.

mtx_cell <- readRDS("/suffolk/WorkGenomicsD/ec474/mitox_azimuth/mtx_singlets_mgatk_cellcycleregressed_annotated_wnn_azimuth.rds")
DefaultAssay(mtx_cell) <- "RNA"

table(mtx_cell@meta.data$predicted.celltype.l1.azimuth, mtx_cell@meta.data$predicted.celltype.l1.wnn)
#             B CD4 T CD8 T  DC Mono  NK other T
#   B        71     0     0   0    0   0       0
#   CD4 T     3   278     4   0    0   0       0
#   CD8 T     0    19    79   0    0  10       1
#   DC        0     0     0   3    0   0       0
#   Mono      0     0     0   0    2   0       0
#   NK        0     0     1   0    0 123       0
#   other T   0     7     2   0    0   0      31

# Normalize the data
ifnb <- NormalizeData(ifnb)

# Find DE features between CD16 Mono and CD1 Mono
Idents(mtx_cell) <- "predicted.celltype.l1.azimuth"
tcell.de.markers <- FindMarkers(mtx_cell, ident.1 = "CD8 T", ident.2 = "CD4 T")

# view results
head(tcell.de.markers)
#              p_val avg_log2FC pct.1 pct.2    p_val_adj
# CCL5  2.157810e-29   2.697702 0.761 0.193 7.897801e-25
# AOAH  3.200153e-23   4.027304 0.431 0.039 1.171288e-18
# NKG7  2.225913e-22   2.697829 0.651 0.186 8.147065e-18
# ZEB2  3.726770e-22   3.972473 0.422 0.042 1.364035e-17
# GNLY  5.790474e-21   2.770826 0.780 0.435 2.119372e-16
# KLRD1 1.590300e-18   4.171805 0.339 0.028 5.820657e-14

#Plot DE genes 
FeaturePlot(mtx_cell, features = "CCL5", reduction="wnn.umap")
ggsave(filename="Tcell_CD8_CD4_CCL5_azimuth.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_Tcell/")


#### WNN
Idents(mtx_cell) <- "predicted.celltype.l1.wnn"
tcell.de.markers <- FindMarkers(mtx_cell, ident.1 = "CD8 T", ident.2 = "CD4 T")

# view results
head(tcell.de.markers)
#             p_val avg_log2FC pct.1 pct.2    p_val_adj
# CCL5  1.934688e-20   2.368291 0.733 0.247 7.081152e-16
# AOAH  3.367857e-18   3.345032 0.407 0.053 1.232669e-13
# ZEB2  6.125417e-17   3.309787 0.407 0.062 2.241964e-12
# NKG7  2.091262e-16   2.241562 0.640 0.220 7.654228e-12
# GNLY  2.515744e-14   2.577618 0.733 0.464 9.207874e-10
# KLRD1 9.163231e-13   3.441856 0.291 0.039 3.353834e-08

#OUTPUT
    # p_val : p-value (unadjusted)
    # avg_log2FC : log fold-change of the average expression between the two groups. 
    #Positive values indicate that the feature is more highly expressed in the first group.
    # pct.1 : The percentage of cells where the feature is detected in the first group
    # pct.2 : The percentage of cells where the feature is detected in the second group
    # p_val_adj : Adjusted p-value, based on Bonferroni correction using all features in the dataset.


#Plot DE genes 
FeaturePlot(mtx_cell, features = "CCL5", reduction="wnn.umap")
ggsave(filename="Tcell_CD8_CD4_CCL5_wnn.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_Tcell/")

### DOT PLOT
marker_gene_list <- c("CD19", "CD4", "CD8A", "FOXP3", "CD86", "CD14")

# B Cells
#     CD19: A pan-B cell marker involved in B cell development, differentiation, and signaling.
#     CD20 (MS4A1): Another pan-B cell marker that is commonly used for the identification of B cells in flow cytometry and histology.
#     CD79A: Part of the B cell receptor complex, essential for B cell activation and signaling.

# T CD4+ Cells
#     CD4: The primary marker for T helper cells, involved in immune signaling and interaction with MHC class II molecules.
#     FOXP3: Marker for regulatory T cells (a subset of CD4+ cells), crucial for maintaining immune tolerance.
#     IL2RA (CD25): High-affinity receptor for IL-2, often expressed on activated T cells and regulatory T cells.

# CD8+ Cells
#     CD8A: Core marker for cytotoxic T cells, involved in interaction with MHC class I molecules.
#     GZMB (Granzyme B): A key effector molecule in the cytotoxic response, involved in inducing apoptosis in target cells.
#     PRF1 (Perforin): Works in concert with granzymes to mediate cell killing by creating pores in the target cell membrane.

# Dendritic Cells
#     CD11c (ITGAX): Commonly used marker for myeloid dendritic cells, involved in antigen presentation and immune activation.
#     HLA-DR: MHC class II molecule, crucial for antigen presentation to T cells.
#     CD86: Co-stimulatory molecule necessary for T cell activation.

# Monocytes
#     CD14: A co-receptor for the detection of bacterial lipopolysaccharide (LPS), widely used as a marker for monocytes.
#     CD16 (FCGR3A): Fc receptor often used to differentiate between classical (CD14++CD16-) and non-classical/intermediate (CD14+CD16+) monocytes.
#     CD68: A marker for monocytes/macrophages involved in phagocytosis.

# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
Idents(mtx_cell) <- "predicted.celltype.l1.wnn"
DotPlot(mtx_cell, features = marker_gene_list, cols=c("lightgrey", "blue", "pink", "orange", "green", "purple", "yellow")) + RotatedAxis()
ggsave(filename="DotPlot_cells_by_marker_genes.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_Tcell/", bg="white")


#Remove cells that are not 
#Create additional metadata column to compare WNN L1 to Azimuth L1 and if they match, label it "match", and then subset the ones that don't match
mtx_cell@meta.data$cell_annotation_concordant_l1 <- ifelse(mtx_cell@meta.data$predicted.celltype.l1.wnn == mtx_cell@meta.data$predicted.celltype.l1.azimuth, "match", "not")

table(mtx_cell@meta.data$cell_annotation_concordant_l1)
# match   not
#   587    47

Idents(mtx_cell) <- "cell_annotation_concordant_l1"
mtx_annotated <- subset(mtx_cell, idents="match")
# An object of class Seurat
# 279441 features across 587 samples within 8 assays
# Active assay: RNA (36601 features, 2000 variable features)
#  3 layers present: counts, data, scale.data
#  7 other assays present: ATAC, HTO, mito, SCT, prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
#  14 dimensional reductions calculated: pca.genes, pca.cellcycle, pca.aftercell, umap.aftercell, umap.aftercell_05, umap.aftercell_07, 
#  umap.aftercell_15, umap.aftercell_2, lsi, umap.atac, wnn.umap, ref.spca, ref.umap, integrated_dr

saveRDS(mtx_annotated, file="/suffolk/WorkGenomicsD/ec474/mitox_azimuth/mtx_singlets_mgatk_cellcycleregressed_annotated_concordant.rds")


##FIND MARKER GENES ON CONCORDANT ASSAY 
#"predicted.celltype.l1"

table(mtx_cell@meta.data$predicted.celltype.l1)

    #   B   CD4 T   CD8 T      DC    Mono      NK other T
    #  71     278      79       3       2     123      31

#Create script that compares one cell type against all others in a seurat object and saves the top 10 DEG
#DEG cannot be performed if fewer than 3 cells in a category

cell.de.markers <- function(seurat_object) {
    DefaultAssay(seurat_object) <- "RNA"
    Idents(seurat_object) <- seurat_object@meta.data$predicted.celltype.l1
    labels <- unique(seurat_object@meta.data$predicted.celltype.l1)
    results <- list()

    for (cell_type in labels) {
    # Check the number of cells in the current cell type
        cells_in_group <- WhichCells(seurat_object, idents = cell_type)
        if (length(cells_in_group) < 3) {
            message(paste("Skipping cell type", cell_type, "due to fewer than 3 cells"))
            next
        }
        if (length(cells_in_group) >= 3) {
    # Perform differential expression analysis
                f <- FindMarkers(seurat_object, ident.1= cell_type, ident.2=NULL)
                f_10 <- as.data.frame(f) %>% slice_head(n=10)
                file_path <- paste0("/suffolk/WorkGenomicsD/ec474/mitox_deg/mitox_", cell_type, "_top10_deg", ".csv")
                write.csv(f_10, file=file_path)
                results[[cell_type]] <- f_10}
                }
    return (results)
}

marker_gene_list <- c("CD19", "CD4", "CD8A", "FOXP3", "CD86", "CD14")
marker_gene_list_2 <- c("INPP4B", "SLC4A10", "LINC02446", "GNLY", "BANK1", "TNFRSF21")

cell.de.markers(mtx_cell)
# Skipping cell type Mono due to fewer than 3 cells
# $`CD4 T`
#               p_val avg_log2FC pct.1 pct.2    p_val_adj
# INPP4B 8.178676e-45  2.0911419 0.878 0.314 2.993477e-40
# NKG7   1.496549e-41 -3.5108461 0.191 0.693 5.477518e-37
# ZEB2   1.966270e-37 -4.6434550 0.043 0.521 7.196746e-33
# GNLY   9.995466e-37 -3.8807786 0.435 0.761 3.658441e-32
# LYN    4.844212e-34 -3.9555214 0.036 0.492 1.773030e-29
# CAMK4  8.201539e-30  1.5505939 0.806 0.353 3.001845e-25
# LEF1   1.596491e-29  2.1732940 0.662 0.233 5.843318e-25
# IL7R   6.560584e-28  1.3883464 0.835 0.388 2.401239e-23
# EEF1A1 2.116855e-27  0.7922189 0.996 0.981 7.747901e-23
# CCL5   5.941735e-27 -2.5991200 0.198 0.599 2.174734e-22

# $`other T`
#                p_val avg_log2FC pct.1 pct.2    p_val_adj
# SLC4A10 3.170237e-26   4.229313 0.516 0.040 1.160338e-21
# KLRG1   5.308301e-22   3.391446 0.645 0.088 1.942891e-17
# IL23R   2.042258e-20   6.049274 0.226 0.005 7.474869e-16
# A2M     5.289438e-14   2.983913 0.710 0.196 1.935987e-09
# GZMK    7.942554e-14   2.444524 0.581 0.108 2.907054e-09
# PLD1    2.737231e-13   4.581955 0.226 0.014 1.001854e-08
# ADAM12  3.934080e-13   4.607893 0.161 0.005 1.439913e-08
# BHLHE40 8.634272e-11   2.653490 0.452 0.088 3.160230e-06
# LTK     1.595954e-10   4.434685 0.161 0.009 5.841352e-06
# TRGC2   2.074380e-10   3.473917 0.355 0.058 7.592437e-06

# $`CD8 T`
#                   p_val avg_log2FC pct.1 pct.2    p_val_adj
# LINC02446  1.710238e-15  4.0219694 0.228 0.022 6.259641e-11
# CD8B       3.463784e-13  2.7849163 0.266 0.041 1.267780e-08
# CCL5       6.274452e-13  1.1919859 0.772 0.352 2.296512e-08
# SGCD       5.047920e-11  3.0955643 0.203 0.028 1.847589e-06
# CD8A       1.342400e-09  2.8982706 0.228 0.043 4.913320e-05
# NIBAN1     1.734484e-07  1.7902451 0.418 0.187 6.348385e-03
# LAG3       2.931278e-07  3.3318779 0.101 0.010 1.072877e-02
# HERC5      1.383922e-06  1.7579405 0.215 0.057 5.065294e-02
# AL360270.1 2.196052e-06  4.2731603 0.076 0.006 8.037769e-02
# PDE3B      3.228848e-06  0.9304706 0.810 0.604 1.181791e-01

# $NK
#               p_val avg_log2FC pct.1 pct.2    p_val_adj
# GNLY   3.754000e-66   3.759330 1.000 0.502 1.374002e-61
# NKG7   6.608611e-63   3.273196 0.967 0.319 2.418818e-58
# PRF1   2.045341e-58   4.347986 0.715 0.075 7.486152e-54
# GZMB   3.046644e-56   4.118916 0.732 0.097 1.115102e-51
# SPON2  7.303129e-52   4.376302 0.610 0.050 2.673018e-47
# KLRF1  4.544938e-46   4.702326 0.537 0.041 1.663493e-41
# TYROBP 1.170723e-44   3.288043 0.667 0.097 4.284963e-40
# GZMA   3.264301e-42   2.751177 0.789 0.205 1.194767e-37
# ZEB2   3.537117e-41   2.697822 0.756 0.172 1.294620e-36
# KLRD1  1.344550e-39   3.455727 0.593 0.084 4.921188e-35

# $B
#                  p_val avg_log2FC pct.1 pct.2    p_val_adj
# BANK1    6.374814e-104   6.894141 0.958 0.025 2.333246e-99
# ARHGAP24  7.918673e-92   6.104387 0.845 0.017 2.898313e-87
# ADAM28    1.528259e-76   5.732369 0.761 0.023 5.593580e-72
# COL19A1   1.921932e-73   7.927664 0.634 0.006 7.034464e-69
# COBLL1    6.271114e-71   5.862247 0.634 0.008 2.295290e-66
# AFF3      5.094873e-70   4.702259 0.972 0.118 1.864774e-65
# FCRL1     2.333991e-68   5.837433 0.648 0.014 8.542639e-64
# EBF1      7.067392e-67   7.388732 0.606 0.010 2.586736e-62
# IGHM      8.204950e-66   5.402026 0.803 0.060 3.003094e-61
# MS4A1     6.107094e-64   5.910961 0.634 0.017 2.235258e-59

# $DC
#                    p_val avg_log2FC pct.1 pct.2     p_val_adj
# TNFRSF21   2.602716e-129  11.500992 1.000 0.000 9.526200e-125
# FLT3        1.157957e-86  11.768443 0.667 0.000  4.238239e-82
# LILRA4      1.157957e-86  11.034024 0.667 0.000  4.238239e-82
# NRP1        2.149473e-58   9.275621 0.667 0.002  7.867285e-54
# DACH1       5.246202e-58   7.135989 0.667 0.002  1.920163e-53
# TUBB6       9.526109e-57   8.352714 1.000 0.007  3.486651e-52
# LINC02818   4.257350e-44  10.306096 0.333 0.000  1.558233e-39
# AL391832.2  4.257350e-44  10.002126 0.333 0.000  1.558233e-39
# AL451007.3  4.257350e-44  10.002126 0.333 0.000  1.558233e-39
# AL451007.2  4.257350e-44  10.002126 0.333 0.000  1.558233e-39

#Dot plot of DEG
marker_gene_list_2 <- c("INPP4B", "SLC4A10", "LINC02446", "GNLY", "BANK1", "TNFRSF21")
Idents(mtx_cell) <- "predicted.celltype.l1.wnn"
DotPlot(mtx_cell, features = marker_gene_list_2, cols=c("lightgrey", "blue", "pink", "orange", "green", "purple", "yellow")) + RotatedAxis()
ggsave(filename="DotPlot_cells_by_marker_genes_in_deg.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_Tcell/", bg="white")


#If the ident.2 parameter is omitted or set to NULL, FindMarkers() will test for differentially expressed features between the group specified by ident.1 and all other cells. 
# Find differentially expressed features between CD14+ Monocytes and all other cells, only
# search for positive markers
monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = NULL, only.pos = TRUE)

# view results
head(monocyte.de.markers)

#Perform DE analysis within the same cell type across conditions
#ince this dataset contains treatment information (control versus stimulated with interferon-beta), we can also ask what genes change in different conditions for cells of the same type. 
#First, we create a column in the meta.data slot to hold both the cell type and treatment information and switch the current Idents to that column. 
#Then we use FindMarkers() to find the genes that are different between control and stimulated CD14 monocytes.

ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"

mono.de <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
head(mono.de, n = 10)

#However, the p-values obtained from this analysis should be interpreted with caution, 
#because these tests treat each cell as an independent replicate and 
#ignore inherent correlations between cells originating from the same sample. 
#Such analyses have been shown to find a large number of false positive associations, 
#as has been demonstrated by Squair et al., 2021, Zimmerman et al., 2021, Junttila et al., 2022, and others. Below, we show how pseudobulking can be used to account for such within-sample correlation.
