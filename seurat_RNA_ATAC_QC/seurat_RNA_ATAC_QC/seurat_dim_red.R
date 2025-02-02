#Variable genes identification, dimensionality reduction with PCA and clustering
library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(Signac)
library(sctransform)
library(clustree)

mtx_qc <- readRDS("/Users/ec474/Mount/Suffolk/RawGenomicsF/ec474/mitox_sce/mtx_qc.rds")

#An object of class Seurat 
#134807 features across 9907 samples within 2 assays 
#Active assay: RNA (36601 features, 0 variable features)
#1 other assay present: ATAC
##1 layer present: counts

head(mtx_qc@meta.data)

# RNA analysis 
DefaultAssay(mtx_qc) <- "RNA"

#Normalise
mtx_norm <- NormalizeData(mtx_qc)

#Find Variable Features
mtx_var <- FindVariableFeatures(mtx_norm, selection.method = 'vst', assay="RNA", layer="count")
top10 <- head(VariableFeatures(mtx_var), 10)

plot1 <- VariableFeaturePlot(mtx_var)
plot2 <- LabelPoints(plot=plot1, points = top10, repel = TRUE)
plot1+plot2

#Scale
mtx_scaled <- ScaleData(mtx_var, features = rownames(mtx_var))

#PCA
mtx_scaled <- RunPCA(mtx_scaled, features = VariableFeatures(mtx_scaled))
print(mtx_scaled[['pca']], dims=1:5, nfeatures=5)

#Visualise PCA
DimPlot(mtx_scaled, reduction='pca') + NoLegend()

VizDimLoadings(mtx_scaled, dims= 1:2, reduction='pca')

DimHeatmap(mtx_scaled, dims=1, cells=500, balanced=TRUE)

#Determine the dimensionality of the data
ElbowPlot(mtx_scaled, ndims = 50)

#Find Neighbours and cluster
mtx_cluster <- FindNeighbors(mtx_scaled, dims= 1:40)

#Intermediate step with clustree to find optimal resolution (do later)
#Cluster
mtx_cluster <- FindClusters(mtx_cluster, resolution = 0.5)

#UMAP, before cell cycle scoring 
mtx_cluster <- RunUMAP(mtx_cluster, dims = 1:40, verbose = FALSE)
DimPlot(mtx_cluster, reduction='umap', label=TRUE)

#Cell cycle scoring ####
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

mtx_scaled_cellcycle <- CellCycleScoring(mtx_cluster, s.features = s.genes, g2m.features = g2m.genes, new.assay.name = 'NORM_CELL', set.ident = TRUE)

head(mtx_scaled_cellcycle@meta.data)

#Visualise distribution of cell cycle markers across features
RidgePlot(mtx_scaled_cellcycle, features = c("JCHAIN", "GNLY", "TOP2A", "MCM6", "MKI67"))

#Plts cell cycles 

# SCTransform (Normalisation, Variable Features, Scaling) regressing out cell cycle stage
mtx_qc <- SCTransform(mtx_qc, verbose = FALSE, vars.to.regress = ') 



#RunPCA()

#ClusterTree to help you find hyperparameter n for KNN
#install.packages("clustree")


#Neighbours
#Clusters

#UMAP
#RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#Cell Annotation
