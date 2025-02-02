
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(sctransform)
library(clustree)
library(patchwork)

mtx_singlet <- readRDS("/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlet_HASH3_RNA_ATAC_QC.rds")

vln <- VlnPlot(mtx_singlet, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent.mt"), ncol = 5, log = TRUE, pt.size = 0.01) 
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "Demultiplexed sample QC")
ggsave("VlnPlot_singlet_points.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_qc_200cells/", width=30, units="cm")

vln <- VlnPlot(mtx_singlet, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent.mt"), ncol = 5, log = TRUE, pt.size = 0) 
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "Demultiplexed sample QC")
ggsave("VlnPlot_singlet_clean.pdf", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_qc_200cells/", width=30, units="cm")
ggsave("VlnPlot_singlet_clean.jpg", path="/suffolk/Homes/ec474/scripts/seurat_HTODemux/figures_qc_200cells/", width=30, units="cm")