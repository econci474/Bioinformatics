library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(patchwork)

mtx <- readRDS("/suffolk/RawGenomicsF/ec474/mitox_sce/mtx.rds")

#mtx.original
#An object of class Seurat 
#36601 features across 10093 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
#1 layer present: counts

#mtx with ATAC
#An object of class Seurat 

#134807 features across 10093 samples within 2 assays 
#Active assay: RNA (36601 features, 0 variable features)
#1 layer present: counts
#1 other assay present: ATAC

#Feature Scatter plots to visualise relationship btw e.g. nCount and nFeature
plot1 <- FeatureScatter(mtx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(mtx, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")
plot1+plot2
ggsave("FeatureScatter_plots.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/", width=30, units = "cm")
ggsave("FeatureScatter_plots.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/", width=30, units = "cm")

#Add feature for percentage mitochondrial DNA 
head(mtx@meta.data)
mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern="^MT-")

#Violin plot of number of detected molecules for each modality as well as mitochondrial percentage
vln <- VlnPlot(mtx, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) 
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "Before QC")
ggsave("VlnPlot_before.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# nCount_ATAC, Lose 37 cells
VlnPlot(mtx, assay="RNA", features="nCount_ATAC", pt.size = 0) + NoLegend() + geom_hline(yintercept=150000, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_before_nCount_ATAC.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_before_nCount_ATAC.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# nFeature_ATAC, Lose 23 cells
VlnPlot(mtx, features="nFeature_ATAC", pt.size = 0) + NoLegend() + geom_hline(yintercept=35000, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_before_nFeature_ATAC.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_before_nFeature_ATAC.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# nCount_RNA, Lose 53 cells from top end, Lose 68 from bottom end with threshold of 100
VlnPlot(mtx, features="nCount_RNA", pt.size = 0) + NoLegend() + geom_hline(yintercept=25000, color="red") + geom_hline(yintercept=100, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_before_nCount_RNA.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_before_nCount_RNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# nFeature_RNA, Lose 21 cells from top end and 121 from bottom end with a threshold of 280
VlnPlot(mtx, features="nFeature_RNA", pt.size = 0) + NoLegend() + geom_hline(yintercept=7500, color="red") + geom_hline(yintercept=280, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_before_nFeature_RNA.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_before_nFeature_RNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# percent.mt, Lose 32 cells
VlnPlot(mtx, features="percent.mt", pt.size=0) + NoLegend() + geom_hline(yintercept=5, colour="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_before_percent_mt.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_before_percent_mt.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# Total estimated cells lost to QC (assuming non-overlapping) = 37 + 23 + 53 + 68 + 21 + 121 + 32 cells 
# which is equal to 323 cells out of 10093 cells = 3.2% poor quality cells

mtx_qc <- subset(
  x = mtx,
  subset = 
    nCount_ATAC < 150000 & 
    nFeature_ATAC < 35000 &
    nCount_RNA < 25000 &
    nCount_RNA > 100 &
    nFeature_RNA < 7500 &
    nFeature_RNA > 280 &
    percent.mt < 5)

#mtx after QC
#An object of class Seurat 
#134807 features across 9907 samples within 2 assays 
#Active assay: RNA (36601 features, 0 variable features)
#1 layer present: counts
#1 other assay present: ATAC

# Total cells lost = 10093 - 9907 = 186, which is 1.8% of all cells  

#Violin plot of number of detected molecules for each modality as well as mitochondrial percentage
vln <- VlnPlot(mtx_qc, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend() + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "After QC")
ggsave("VlnPlot_after.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# Not cutting bottom end of ATAC. What is the significance of zero open chromatin regions? 
# nCount_ATAC, Lose 37 cells
VlnPlot(mtx_qc, features="nCount_ATAC", pt.size = 0) + NoLegend() + geom_hline(yintercept=150000, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_after_nCount_ATAC.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_after_nCount_ATAC.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# nFeature_ATAC, Lose 23 cells
VlnPlot(mtx_qc, features="nFeature_ATAC", pt.size = 0) + NoLegend() + geom_hline(yintercept=35000, colour="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_after_nFeature_ATAC.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_after_nFeature_ATAC.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# nCount_RNA, Lose 53 cells from top end, Lose 68 from bottom end with threshold of 100
VlnPlot(mtx_qc, features="nCount_RNA", pt.size = 0) + NoLegend() + geom_hline(yintercept=25000, color="red") + geom_hline(yintercept=50, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("VlnPlot_after_nCount_RNA.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_after_nCount_RNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# nFeature_RNA, Lose 21 cells from top end and 121 from bottom end with a threshold of 280
VlnPlot(mtx_qc, features="nFeature_RNA", pt.size = 0) + NoLegend() + geom_hline(yintercept=7500, color="red") + geom_hline(yintercept=280, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_after_nFeature_RNA.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_after_nFeature_RNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

# percent.mt, Lose 32 cells
VlnPlot(mtx_qc, features="percent.mt", pt.size=0) + NoLegend() + geom_hline(yintercept=5, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_after_percent_mt.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave("VlnPlot_after_percent_mt.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")


# Feature plot
DefaultAssay(mtx_qc) <- "RNA"
mtx_qc <- NormalizeData(mtx_qc)
mtx_qc <- FindVariableFeatures(mtx_qc, selection.method = "vst")
mtx_qc <- ScaleData(mtx_qc, features = rownames(mtx_qc))
mtx_qc <- RunPCA(mtx_qc, features = VariableFeatures(object=mtx_qc), reduction.name="pca_qc", reduction.key="qcPCA_")

ElbowPlot(mtx_qc, reduction="pca_qc", ndims=50)
ggsave("Elbowplot_after.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

mtx_qc <- FindNeighbors(mtx_qc, reduction="pca_qc", dims=1:50)

#Best clustering resolution
# Define the resolution values and loop over them
R_values <- c(0.2, 0.5, 0.7, 1, 1.5, 2)

DefaultAssay(mtx_qc) <- "RNA"
for (R in R_values) {
  mtx_qc <- FindClusters(mtx_qc, verbose = FALSE, resolution = R, random.seed=1)
  R_key <- gsub("\\.", "", as.character(R))
  mtx_qc <- RunUMAP(mtx_qc, dims = 1:50, verbose = FALSE, reduction = "pca_qc",
                      reduction.name=paste0("umap.clus_", R_key), reduction.key = paste0("clusUMAP", R_key, "_"))
  plot <- DimPlot(mtx_qc, reduction = paste0("umap.clus_", R_key), label = TRUE) + 
    labs(title = paste("UMAP with clustering resolution of", R))
  ggsave(plot, filename = paste0("umap_after_cluster_", R_key,".pdf"), path = "/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_clustering/")
}

#find prefix 
head(mtx_qc@meta.data)

#Visualise and save clustree
library(clustree)
clustree <- clustree(mtx_qc, prefix="RNA_snn_res.")
ggsave(clustree, filename="clustree.pdf", path= "/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_clustering/", height=10)
ggsave(clustree, filename="clustree.jpg", path= "/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_clustering/", height=10)

#Best clustering is of resolution = 1
mtx_qc <- FindNeighbors(mtx_qc, reduction="pca_qc", dims=1:50)
mtx_qc <- FindClusters(mtx_qc, verbose = FALSE, resolution = 1, random.seed=1)
plot <- DimPlot(mtx_qc, reduction="umap.clus_1", label=TRUE) + labs(title = paste("UMAP with clustering resolution of 1"))
ggsave(plot, filename="UMAP_after.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")
ggsave(plot, filename="UMAP_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures/")

#Featureplot on best clustering
#nCount_ATAC
plot <- FeaturePlot(mtx_qc, feature="nCount_ATAC", reduction="umap.clus_1") + labs(title = paste("FeaturePlot of nCount_ATAC"))
ggsave(plot, filename="UMAP_after_nCount_ATAC.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")
ggsave(plot, filename="UMAP_after_nCount_ATAC.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")

#nFeature_ATAC
plot <- FeaturePlot(mtx_qc, feature="nFeature_ATAC", reduction="umap.clus_1") + labs(title = paste("FeaturePlot of nFeature_ATAC"))
ggsave(plot, filename="UMAP_after_nFeature_ATAC.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")
ggsave(plot, filename="UMAP_after_nFeature_ATAC.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")

#nCount_RNA
plot <- FeaturePlot(mtx_qc, feature="nCount_RNA", reduction="umap.clus_1") + labs(title = paste("FeaturePlot of nCount_RNA"))
ggsave(plot, filename="UMAP_after_nCount_RNA.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")
ggsave(plot, filename="UMAP_after_nCount_RNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")

#nFeature_RNA
plot <- FeaturePlot(mtx_qc, feature="nFeature_RNA", reduction="umap.clus_1") + labs(title = paste("FeaturePlot of nFeature_RNA"))
ggsave(plot, filename="UMAP_after_nFeature_RNA.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")
ggsave(plot, filename="UMAP_after_nFeature_RNA.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")

#percent_mt
plot <- FeaturePlot(mtx_qc, feature="percent.mt", reduction="umap.clus_1") + labs(title = paste("FeaturePlot of percent_mt"))
ggsave(plot, filename="UMAP_after_percent_mt.pdf", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")
ggsave(plot, filename="UMAP_after_percent_mt.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_featureplots/")


saveRDS(object=mtx_qc, file="/suffolk/RawGenomicsF/ec474/mitox_sce/mtx_qc.rds")
