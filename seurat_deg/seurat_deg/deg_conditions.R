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

mtx_cell <- readRDS("/suffolk/WorkGenomicsD/ec474/mitox_azimuth/mtx_singlets_mgatk_cellcycleregressed_annotated_concordant.rds")
DefaultAssay(mtx_cell) <- "RNA"

# create a new column to annotate sample-condition-celltype in the single-cell dataset
mtx_cell@meta.data$hash.celltype <- paste0(mtx_cell@meta.data$HTO_maxID, "-", mtx_cell@meta.data$predicted.celltype.l1)

table(mtx_cell@meta.data$hash.celltype)
# NPC-HASH1-B   NPC-HASH1-CD4 T   NPC-HASH1-CD8 T      NPC-HASH1-DC
#                21               104                21                 2
#    NPC-HASH1-Mono      NPC-HASH1-NK NPC-HASH1-other T       NPC-HASH2-B
#                 1                41                14                18
#   NPC-HASH2-CD4 T   NPC-HASH2-CD8 T      NPC-HASH2-DC      NPC-HASH2-NK
#                53                25                 1                39
# NPC-HASH2-other T       NPC-HASH3-B   NPC-HASH3-CD4 T   NPC-HASH3-CD8 T
#                 8                32               121                33
#    NPC-HASH3-Mono      NPC-HASH3-NK NPC-HASH3-other T
#                 1                43                 9

# generate violin plot 
Idents(mtx_cell) <- "hash.celltype"

#DE between same cell type across HASH 1 (control) and HASH 2 (before intervention) 
cell.de.markers_hash1_hash2 <- function(seurat_object) {
    DefaultAssay(seurat_object) <- "RNA"
    Idents(seurat_object) <- seurat_object@meta.data$hash.celltype
    labels <- unique(seurat_object@meta.data$hash.celltype)
    results <- list()
    base_cell_types <- unique(seurat_object@meta.data$predicted.celltype.l1)

    for (base_cell_type in base_cell_types) {
        cell_type_hash1 <- paste0("NPC-HASH1-", base_cell_type)
        cell_type_hash2 <- paste0("NPC-HASH2-", base_cell_type)
        
        # Check if both cell types exist in the data
        if (cell_type_hash1 %in% labels & cell_type_hash2 %in% labels) {
            # Check the number of cells in each group
            cells_in_group1 <- WhichCells(seurat_object, idents = cell_type_hash1)
            cells_in_group2 <- WhichCells(seurat_object, idents = cell_type_hash2)
            
            if (length(cells_in_group1) >= 3 & length(cells_in_group2) >= 3) {
                
                # Perform differential expression analysis
                f <- FindMarkers(seurat_object, ident.1 = cell_type_hash1, ident.2 = cell_type_hash2)
                f_10 <- as.data.frame(f) %>% slice_head(n = 10)
                file_path <- paste0("/suffolk/WorkGenomicsD/ec474/mitox_deg/mitox_", base_cell_type, "_HASH1_vs_HASH2_top10_deg.csv")
                write.csv(f_10, file = file_path)
                results[[paste0(base_cell_type, "_HASH1_vs_HASH2")]] <- f_10
            } else {
                message(paste("Skipping comparison for cell type", base_cell_type, "due to fewer than 3 cells in one or both groups"))
            }
        } else {
            message(paste("Skipping comparison for cell type", base_cell_type, "because one of the groups does not exist in the data"))
        }
    }
    
    return(results)   
}

cell.de.markers_hash1_hash2(mtx_cell)
# Skipping comparison for cell type Mono because one of the groups does not exist in the data
# Skipping comparison for cell type DC due to fewer than 3 cells in one or both groups
# $`CD4 T_HASH1_vs_HASH2`
#                   p_val avg_log2FC pct.1 pct.2 p_val_adj
# FOS        1.387849e-05 -1.0170266 0.317 0.698 0.5079667
# TNFAIP3    1.702548e-05 -1.3711254 0.288 0.660 0.6231497
# TRAF3      4.836505e-05 -1.1077222 0.135 0.453 1.0000000
# YY1        9.320176e-05 -1.5992157 0.125 0.396 1.0000000
# SCFD1      1.113547e-04 -1.3410184 0.212 0.509 1.0000000
# HNRNPDL    1.211814e-04 -0.9610802 0.413 0.698 1.0000000
# AL035634.1 1.612223e-04 -6.0816975 0.000 0.132 1.0000000
# RPL41      2.261411e-04  0.3460652 1.000 1.000 1.0000000
# NSMCE2     2.865583e-04 -1.2040219 0.221 0.528 1.0000000
# RORA-AS1   3.090786e-04 -1.4954575 0.202 0.472 1.0000000

# $`other T_HASH1_vs_HASH2`
#                p_val avg_log2FC pct.1 pct.2 p_val_adj
# EIF1B    0.001343150  -6.083054 0.000 0.625         1
# ESYT2    0.001343150  -5.526745 0.000 0.625         1
# SLC4A10  0.003297157   3.668646 0.786 0.125         1
# ITK      0.005284594  -4.905501 0.000 0.500         1
# HIGD2A   0.005284594  -5.545949 0.000 0.500         1
# GSTK1    0.005284594  -5.728096 0.000 0.500         1
# VCP      0.005284594  -5.652965 0.000 0.500         1
# RHOG     0.005284594  -5.752674 0.000 0.500         1
# ATP5F1B  0.005284594  -5.023988 0.000 0.500         1
# CALCOCO2 0.005284594  -5.366139 0.000 0.500         1

# $`CD8 T_HASH1_vs_HASH2`
#                p_val avg_log2FC pct.1 pct.2 p_val_adj
# FAM172A 0.0002620282 -1.7988920 0.286  0.88         1
# FOSB    0.0003523769 -3.8748233 0.048  0.56         1
# PNPLA8  0.0006990189 -5.8146181 0.000  0.44         1
# PDCD4   0.0010731827 -1.5966895 0.333  0.80         1
# BIN2    0.0011324448 -2.7048942 0.095  0.56         1
# EIF3H   0.0011324448  2.8923666 0.571  0.16         1
# TMSB4X  0.0014953172  0.7575151 1.000  1.00         1
# RNMT    0.0016333461 -2.8566187 0.095  0.52         1
# PIK3CA  0.0020585112 -2.8196174 0.048  0.48         1
# FCMR    0.0021442337  5.7337822 0.333  0.00         1

# $NK_HASH1_vs_HASH2
#                p_val avg_log2FC pct.1 pct.2 p_val_adj
# HLA-F   0.0002922906   2.167282 0.537 0.154         1
# SPN     0.0004992105  -2.132214 0.073 0.436         1
# GATAD2B 0.0005802534   5.884261 0.268 0.000         1
# STYX    0.0012454211  -5.614359 0.000 0.231         1
# LUC7L   0.0012593088  -1.565827 0.073 0.410         1
# TSPO    0.0014947594  -2.663979 0.073 0.359         1
# RIPOR2  0.0016982661   1.250230 0.756 0.538         1
# TBCK    0.0018077262  -1.615540 0.073 0.385         1
# ARMH4   0.0021379814   6.010144 0.220 0.000         1
# MOB1A   0.0021587620   3.279280 0.317 0.051         1

# $B_HASH1_vs_HASH2
#                  p_val avg_log2FC pct.1 pct.2 p_val_adj
# MAPK8IP3   0.000334538  -5.691894 0.000 0.500         1
# LINC-PINT  0.002116877  -5.171936 0.000 0.389         1
# U2SURP     0.002399202  -2.332676 0.143 0.611         1
# ARFGEF1    0.003040537  -1.311790 0.095 0.556         1
# USP3       0.004202101  -1.567144 0.333 0.778         1
# NCK2       0.004247417  -1.681914 0.190 0.667         1
# TCF7       0.005036215  -3.299726 0.048 0.444         1
# AL512603.2 0.005048733  -5.190511 0.000 0.333         1
# PLEKHA1    0.005048733  -5.980728 0.000 0.333         1
# GOLGA8A    0.005048733  -4.990450 0.000 0.333         1



#DE between same cell type across HASH 2 (before) and HASH 3 (after) 
cell.de.markers_hash2_hash3 <- function(seurat_object) {
    DefaultAssay(seurat_object) <- "RNA"
    Idents(seurat_object) <- seurat_object@meta.data$hash.celltype
    labels <- unique(seurat_object@meta.data$hash.celltype)
    results <- list()
    base_cell_types <- unique(seurat_object@meta.data$predicted.celltype.l1)

    for (base_cell_type in base_cell_types) {
        cell_type_hash2 <- paste0("NPC-HASH2-", base_cell_type)
        cell_type_hash3 <- paste0("NPC-HASH3-", base_cell_type)
        
        # Check if both cell types exist in the data
        if (cell_type_hash2 %in% labels & cell_type_hash3 %in% labels) {
            # Check the number of cells in each group
            cells_in_group2 <- WhichCells(seurat_object, idents = cell_type_hash2)
            cells_in_group3 <- WhichCells(seurat_object, idents = cell_type_hash3)
            
            if (length(cells_in_group2) >= 3 & length(cells_in_group3) >= 3) {
                
                # Perform differential expression analysis
                f <- FindMarkers(seurat_object, ident.1 = cell_type_hash2, ident.2 = cell_type_hash3)
                f_10 <- as.data.frame(f) %>% slice_head(n = 10)
                file_path <- paste0("/suffolk/WorkGenomicsD/ec474/mitox_deg/mitox_", base_cell_type, "_HASH2_vs_HASH3_top10_deg.csv")
                write.csv(f_10, file = file_path)
                results[[paste0(base_cell_type, "_HASH2_vs_HASH3")]] <- f_10
            } else {
                message(paste("Skipping comparison for cell type", base_cell_type, "due to fewer than 3 cells in one or both groups"))
            }
        } else {
            message(paste("Skipping comparison for cell type", base_cell_type, "because one of the groups does not exist in the data"))
        }
    }
    
    return(results)   
}

cell.de.markers_hash2_hash3(mtx_cell)

# $`CD4 T_HASH2_vs_HASH3`
#                 p_val avg_log2FC pct.1 pct.2 p_val_adj
# GNAS     0.0000268904   1.086140 0.698 0.405 0.9842154
# ATOX1    0.0003239836   4.208534 0.132 0.008 1.0000000
# RGS19    0.0004147562   3.237092 0.151 0.017 1.0000000
# VASP     0.0005522668   2.602376 0.170 0.025 1.0000000
# ADIPOR1  0.0005865711   2.086007 0.208 0.041 1.0000000
# ATXN7L3B 0.0005865711   2.393280 0.208 0.041 1.0000000
# CAPNS1   0.0007074890   2.144514 0.245 0.066 1.0000000
# SLC44A2  0.0009294453   1.607909 0.264 0.074 1.0000000
# CWF19L2  0.0009833871  -2.343766 0.094 0.314 1.0000000
# GCLC     0.0012923994   3.716688 0.113 0.008 1.0000000

# $`other T_HASH2_vs_HASH3`
#              p_val avg_log2FC pct.1 pct.2 p_val_adj
# PPIB   0.007887206   2.342958 0.875 0.111         1
# WNK1   0.008620586   5.058909 0.625 0.000         1
# ANXA1  0.010548935  -2.447665 0.375 0.889         1
# ITSN2  0.010690074   2.271320 0.750 0.222         1
# TANK   0.019898757  -5.387855 0.000 0.556         1
# CFLAR  0.019898757  -4.800298 0.000 0.556         1
# SETD5  0.019898757  -5.051710 0.000 0.556         1
# ZNF148 0.019898757  -5.048479 0.000 0.556         1
# RREB1  0.019898757  -5.359234 0.000 0.556         1
# JAK2   0.019898757  -4.762513 0.000 0.556         1

# $`CD8 T_HASH2_vs_HASH3`
#                   p_val avg_log2FC pct.1 pct.2 p_val_adj
# STK17B     0.0003371480   1.477820  0.84 0.364         1
# CNTRL      0.0004852009   3.771402  0.40 0.030         1
# GSTK1      0.0008713272   1.811848  0.52 0.091         1
# PCED1B     0.0009507574  -1.800144  0.20 0.636         1
# BIN2       0.0011308785   1.814655  0.56 0.152         1
# AL360091.3 0.0014175972   5.969661  0.28 0.000         1
# C7orf50    0.0014175972   5.332203  0.28 0.000         1
# MIGA1      0.0023669147   3.511019  0.32 0.030         1
# EIF3G      0.0029797266  -5.455174  0.00 0.303         1
# P4HB       0.0030257942   3.361491  0.32 0.030         1

# $NK_HASH2_vs_HASH3
#                 p_val avg_log2FC pct.1 pct.2 p_val_adj
# HLA-F    0.0001495119 -2.0671455 0.154 0.558         1
# LMTK2    0.0002080995  2.3363734 0.436 0.070         1
# SHPRH    0.0003212522  3.0404324 0.333 0.023         1
# SERPINB1 0.0007955379 -6.2383415 0.000 0.256         1
# XPR1     0.0008008366  1.1794990 0.564 0.163         1
# TERF2IP  0.0012357174  2.1122918 0.385 0.070         1
# KRT10    0.0014201598  1.9678047 0.436 0.116         1
# GATAD2B  0.0014868784 -5.8942330 0.000 0.233         1
# RPS8     0.0016439310  0.5169275 0.974 0.907         1
# DIS3L2   0.0019699446  1.3079999 0.308 0.047         1

# $B_HASH2_vs_HASH3
#                 p_val avg_log2FC pct.1 pct.2 p_val_adj
# HBP1     0.0001581032   2.681487 0.611 0.094         1
# ZRSR2    0.0001876651   6.164907 0.389 0.000         1
# TMEM263  0.0006272345   5.845307 0.333 0.000         1
# PRKDC    0.0008714248   1.337797 0.667 0.188         1
# IL2RG    0.0009020614   1.784969 0.611 0.156         1
# FBXW7    0.0010040290   1.379714 0.778 0.375         1
# MAPK8IP3 0.0010536482   1.748682 0.500 0.062         1
# RNF166   0.0011358810   2.923795 0.444 0.062         1
# LSM10    0.0020274265   5.600200 0.278 0.000         1
# COPS4    0.0020274265   5.572527 0.278 0.000         1

#### COMPARE HASH 3 (after) with HASH 2 (before), it should be the same but opposite direction of log2FC

cell.de.markers_hash3_hash2 <- function(seurat_object) {
    DefaultAssay(seurat_object) <- "RNA"
    Idents(seurat_object) <- seurat_object@meta.data$hash.celltype
    labels <- unique(seurat_object@meta.data$hash.celltype)
    results <- list()
    base_cell_types <- unique(seurat_object@meta.data$predicted.celltype.l1)

    for (base_cell_type in base_cell_types) {
        cell_type_hash2 <- paste0("NPC-HASH2-", base_cell_type)
        cell_type_hash3 <- paste0("NPC-HASH3-", base_cell_type)
        
        # Check if both cell types exist in the data
        if (cell_type_hash2 %in% labels & cell_type_hash3 %in% labels) {
            # Check the number of cells in each group
            cells_in_group2 <- WhichCells(seurat_object, idents = cell_type_hash2)
            cells_in_group3 <- WhichCells(seurat_object, idents = cell_type_hash3)
            
            if (length(cells_in_group2) >= 3 & length(cells_in_group3) >= 3) {
                
                # Perform differential expression analysis
                f <- FindMarkers(seurat_object, ident.1 = cell_type_hash3, ident.2 = cell_type_hash2)
                f_10 <- as.data.frame(f) %>% slice_head(n = 10)
                file_path <- paste0("/suffolk/WorkGenomicsD/ec474/mitox_deg/mitox_", base_cell_type, "_HASH3_vs_HASH2_top10_deg.csv")
                write.csv(f_10, file = file_path)
                results[[paste0(base_cell_type, "_HASH3_vs_HASH2")]] <- f_10
            } else {
                message(paste("Skipping comparison for cell type", base_cell_type, "due to fewer than 3 cells in one or both groups"))
            }
        } else {
            message(paste("Skipping comparison for cell type", base_cell_type, "because one of the groups does not exist in the data"))
        }
    }
    
    return(results)   
}

cell.de.markers_hash3_hash2(mtx_cell)
# $`CD4 T_HASH3_vs_HASH2`
#                 p_val avg_log2FC pct.1 pct.2 p_val_adj
# GNAS     0.0000268904  -1.086140 0.405 0.698 0.9842154
# ATOX1    0.0003239836  -4.208534 0.008 0.132 1.0000000
# RGS19    0.0004147562  -3.237092 0.017 0.151 1.0000000
# VASP     0.0005522668  -2.602376 0.025 0.170 1.0000000
# ADIPOR1  0.0005865711  -2.086007 0.041 0.208 1.0000000
# ATXN7L3B 0.0005865711  -2.393280 0.041 0.208 1.0000000
# CAPNS1   0.0007074890  -2.144514 0.066 0.245 1.0000000
# SLC44A2  0.0009294453  -1.607909 0.074 0.264 1.0000000
# CWF19L2  0.0009833871   2.343766 0.314 0.094 1.0000000
# GCLC     0.0012923994  -3.716688 0.008 0.113 1.0000000

# $`other T_HASH3_vs_HASH2`
#              p_val avg_log2FC pct.1 pct.2 p_val_adj
# PPIB   0.007887206  -2.342958 0.111 0.875         1
# WNK1   0.008620586  -5.058909 0.000 0.625         1
# ANXA1  0.010548935   2.447665 0.889 0.375         1
# ITSN2  0.010690074  -2.271320 0.222 0.750         1
# TANK   0.019898757   5.387855 0.556 0.000         1
# CFLAR  0.019898757   4.800298 0.556 0.000         1
# SETD5  0.019898757   5.051710 0.556 0.000         1
# ZNF148 0.019898757   5.048479 0.556 0.000         1
# RREB1  0.019898757   5.359234 0.556 0.000         1
# JAK2   0.019898757   4.762513 0.556 0.000         1

# $`CD8 T_HASH3_vs_HASH2`
#                   p_val avg_log2FC pct.1 pct.2 p_val_adj
# STK17B     0.0003371480  -1.477820 0.364  0.84         1
# CNTRL      0.0004852009  -3.771402 0.030  0.40         1
# GSTK1      0.0008713272  -1.811848 0.091  0.52         1
# PCED1B     0.0009507574   1.800144 0.636  0.20         1
# BIN2       0.0011308785  -1.814655 0.152  0.56         1
# AL360091.3 0.0014175972  -5.969661 0.000  0.28         1
# C7orf50    0.0014175972  -5.332203 0.000  0.28         1
# MIGA1      0.0023669147  -3.511019 0.030  0.32         1
# EIF3G      0.0029797266   5.455174 0.303  0.00         1
# P4HB       0.0030257942  -3.361491 0.030  0.32         1

# $NK_HASH3_vs_HASH2
#                 p_val avg_log2FC pct.1 pct.2 p_val_adj
# HLA-F    0.0001495119  2.0671455 0.558 0.154         1
# LMTK2    0.0002080995 -2.3363734 0.070 0.436         1
# SHPRH    0.0003212522 -3.0404324 0.023 0.333         1
# SERPINB1 0.0007955379  6.2383415 0.256 0.000         1
# XPR1     0.0008008366 -1.1794990 0.163 0.564         1
# TERF2IP  0.0012357174 -2.1122918 0.070 0.385         1
# KRT10    0.0014201598 -1.9678047 0.116 0.436         1
# GATAD2B  0.0014868784  5.8942330 0.233 0.000         1
# RPS8     0.0016439310 -0.5169275 0.907 0.974         1
# DIS3L2   0.0019699446 -1.3079999 0.047 0.308         1

# $B_HASH3_vs_HASH2
#                 p_val avg_log2FC pct.1 pct.2 p_val_adj
# HBP1     0.0001581032  -2.681487 0.094 0.611         1
# ZRSR2    0.0001876651  -6.164907 0.000 0.389         1
# TMEM263  0.0006272345  -5.845307 0.000 0.333         1
# PRKDC    0.0008714248  -1.337797 0.188 0.667         1
# IL2RG    0.0009020614  -1.784969 0.156 0.611         1
# FBXW7    0.0010040290  -1.379714 0.375 0.778         1
# MAPK8IP3 0.0010536482  -1.748682 0.062 0.500         1
# RNF166   0.0011358810  -2.923795 0.062 0.444         1
# LSM10    0.0020274265  -5.600200 0.000 0.278         1
# COPS4    0.0020274265  -5.572527 0.000 0.278         1



###VISUALISATION 

#Volcano plot
library(ggplot2)

plot_volcano <- function(deg_results, comparison_label) {
    # Create a volcano plot
    ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(alpha = 0.4, aes(color = avg_log2FC > 0)) +
        theme_minimal() +
        labs(title = paste("Volcano plot for", comparison_label),
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value") +
        scale_color_manual(values = c("red", "blue"), 
                           name = "Fold Change",
                           labels = c("Downregulated", "Upregulated")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
}

library(EnhancedVolcano)
plot_enhanced_volcano <- function(deg_results, comparison_label) {
    EnhancedVolcano(deg_results,
                    lab = rownames(deg_results),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = paste("Volcano plot for", comparison_label),
                    xlab = bquote(~Log[2]~ "fold change"),
                    ylab = bquote(~-Log[10]~ "adjusted p-value"),
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    pointSize = 3.0,
                    labSize = 3.0,
                    col = c("grey30", "forestgreen", "royalblue", "red2"),
                    legendPosition = 'right',
                    legendLabSize = 12,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5)
}


# Example usage
results_hash1_hash2 <- cell.de.markers_hash1_hash2(mtx_cell)
comparison_label <- "CD4 T_HASH1_vs_HASH2"
deg_results <- results_hash1_hash2[[comparison_label]]
plot_enhanced_volcano(deg_results, comparison_label)

ggsave(filename="Volcano_CD4_T_HASH1_HAHS2.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_hash1_hash2")


results_hash2_hash3 <- FindMarkers(mtx_cell, ident.1 = "NPC-HASH2-CD4 T", ident.2 = "NPC-HASH3-CD4 T")
comparison_label <- "CD4 T_HASH2_vs_HASH3"

EnhancedVolcano(results_hash2_hash3,
                    lab = rownames(results_hash2_hash3),
                    title = paste("Volcano plot for", comparison_label),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    xlab = bquote(~Log[2]~ "fold change"),                    
                    ylab = bquote(~-Log[10]~ "adjusted p-value"),
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    pointSize = 3.0,
                    labSize = 3.0,
                    col = c("grey30", "forestgreen", "royalblue", "red2"),
                    legendPosition = 'right',
                    legendLabSize = 12,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5)

ggsave(filename="Volcano_CD4_T_HASH2_HAHS3.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_hash2_hash3")


#Heatmap 
library(pheatmap)

plot_heatmap <- function(seurat_object, top_genes, comparison_label) {
    # Extract expression data for top genes
    expression_data <- GetAssayData(seurat_object, slot = "data")[top_genes, ]
    
    # Scale the data
    scaled_data <- t(scale(t(expression_data)))
    
    # Create heatmap
    pheatmap(scaled_data, cluster_rows = TRUE, cluster_cols = TRUE, 
             show_rownames = TRUE, show_colnames = FALSE,
             main = paste("Heatmap for", comparison_label))
}

plot_heatmap <- function(seurat_object, top_genes, comparison_label) {
    # Generate heatmap
    heatmap <- DoHeatmap(seurat_object, features = top_genes, 
                         group.by = "hash.celltype",
                         label = TRUE, 
                         size = 5.5,
                         hjust = 0,
                         vjust = 0,
                         angle = 45,
                         slot = "scale.data", 
                         disp.max = 10) +
        ggtitle(paste("Heatmap for", comparison_label) + 
        theme(axis.title.x = element_blank())
        ) 
    return(heatmap)
}

# Example usage
results_hash1_hash2 <- cell.de.markers_hash1_hash2(mtx_cell)
comparison_label <- "CD4 T_HASH1_vs_HASH2"
deg_results <- results_hash1_hash2[[comparison_label]]

top_genes <- rownames(deg_results)
plot_heatmap(mtx_cell, top_genes, comparison_label)
ggsave(filename="HeatMap_CD4_T_HASH1_HAHS2.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_hash1_hash2")

Idents(mtx_cd4) <- "hash.celltype"
mtx_cd4_hash2_hash3 <- subset(mtx_cd4, idents=c("NPC-HASH2-CD4 T", "NPC-HASH3-CD4 T")

DoHeatmap(mtx_cd4_hash2_hash3, features = rownames(results_hash2_hash3)[1:10], 
                         group.by = "hash.celltype",
                         label = TRUE, 
                         size = 5.5,
                         hjust = 0,
                         vjust = 0,
                         angle = 45,
                         slot = "scale.data", 
                         disp.max = 10) +
        ggtitle(paste("Heatmap for", comparison_label)) + 
        theme(axis.title.x = element_blank())

ggsave(filename="Heatmap_CD4_T_HASH2_HAHS3.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_hash2_hash3")


#Vln plot 
plot_violin <- function(seurat_object, gene, comparison_label) {
    VlnPlot(seurat_object, features = gene, group.by = "hash.celltype") +
        theme_minimal() +
        labs(title = paste("Violin plot for", gene, "in", comparison_label),
             x = "Cell Type",
             y = "Expression Level")
}

# Example usage
results_hash1_hash2 <- cell.de.markers_hash1_hash2(mtx_cell)
comparison_label <- "CD4 T_HASH1_vs_HASH2"
deg_results <- results_hash1_hash2[[comparison_label]]
top_genes <- rownames(deg_results)
gene <- top_genes[1] # Plot the first top gene as an example
plot_violin(mtx_cell, gene, comparison_label)

ggsave(filename="Vln_CD4_T_HASH1_HAHS2.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_hash1_hash2")


# Directory to save plots
output_dir <- "/suffolk/WorkGenomicsD/ec474/mitox_deg/top10_deg_plots_hash1_hash2"

# Loop through all results
for (comparison_label in names(results_hash1_hash2)) {
    deg_results <- results_hash1_hash2[[comparison_label]]
    
    # Volcano Plot
    volcano_plot <- plot_volcano(deg_results, comparison_label)
    ggsave(filename = paste0(output_dir, "/", comparison_label, "_volcano.jpg"), plot = volcano_plot)
    
    # Heatmap
    top_genes <- rownames(deg_results)
    heatmap_plot <- plot_heatmap(mtx_cell, top_genes, comparison_label)
    ggsave(filename = paste0(output_dir, "/", comparison_label, "_heatmap.jpg"), plot = heatmap_plot)
    
    # Violin Plot for top 5 genes
    for (gene in top_genes[1:5]) {
        violin_plot <- plot_violin(mtx_cell, gene, comparison_label)
        ggsave(filename = paste0(output_dir, "/", comparison_label, "_violin_", gene, ".jpg"), plot = violin_plot)
    }
}


ggsave(filename="Volcano_CD4_T_HASH1_HAHS2.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_deg/figures_hash1_hash2")
