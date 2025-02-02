#Code from Signac vignette: https://stuartlab.org/signac/articles/mito 
#Load Seurat object
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)

#Load Seurat object for HASH1 patient - this is not the correct one, needs to process this further first but okay for now as a quick test
mtx <- readRDS(file= "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlet_HASH1_RNA_ATAC_QC.rds")

DefaultAssay(mtx) <- "ATAC"
mtx
# An object of class Seurat
# 134811 features across 221 samples within 3 assays
# Active assay: ATAC (98206 features, 0 variable features)
#  2 layers present: counts, data
#  2 other assays present: RNA, HTO

head(rownames(mtx))
# [1] "chr1-9803-10701"    "chr1-181010-181757" "chr1-778288-779195"
# [4] "chr1-816891-817776" "chr1-819983-820697" "chr1-822687-823558"

head(colnames(mtx))
# [1] "AAACGCGCAGCACGTT-1" "AAACGCGCAGGCAAGC-1" "AACAGGATCGCTAGAT-1"
# [4] "AACCGCTCAACGTGCT-1" "AACTACTCAGGAACCA-1" "AAGAACAGTGTTTCAC-1"

#Load the data from MGATK output: /suffolk/WorkGenomicsD/ec474/mgatk
# load mgatk output
mito.data <- ReadMGATK(dir="/suffolk/WorkGenomicsD/ec474/mgatk/mgatk_hash1/final")
# Reading allele counts
# Reading metadata
# Building matrices

# create an assay
mito <- CreateAssayObject(counts = mito.data$counts)
mito
# Assay data with 132552 features for 221 cells
# First 10 features:
#  A-1-fwd, A-2-fwd, A-3-fwd, A-4-fwd, A-5-fwd, A-6-fwd, A-7-fwd, A-8-fwd,
# A-9-fwd, A-10-fwd
#Letter-position-strand

head(rownames(mito))
# [1] "A-1-fwd" "A-2-fwd" "A-3-fwd" "A-4-fwd" "A-5-fwd" "A-6-fwd"

head(colnames(mito))
# [1] "ATTGGTTCATGAAGTA-1" "CAAACTGGTCAGGCAT-1" "CAAACTGGTCTAACAG-1"
# [4] "CAAGAACCATATAACC-1" "CAAGGGAGTAAGGTTT-1" "CAAGTAACAAACCCTA-1"

# Subset to cell present in the scATAC-seq assay
mito <- subset(mito, cells = Cells(mtx))

# add assay and metadata to the seurat object
mtx[["mito"]] <- mito

mtx
# An object of class Seurat
# 267363 features across 221 samples within 4 assays
# Active assay: ATAC (98206 features, 0 variable features)
#  2 layers present: counts, data
#  3 other assays present: RNA, HTO, mito

mtx <- AddMetaData(mtx, metadata = mito.data$depth[Cells(mito), ], col.name = "mtDNA_depth")

colnames(mtx@meta.data)
#  [1] "orig.ident"                        "nCount_RNA"
#  [3] "nFeature_RNA"                      "nCount_ATAC"
#  [5] "nFeature_ATAC"                     "percent.mt"
#  [7] "nucleosome_signal"                 "nucleosome_percentile"
#  [9] "TSS.enrichment"                    "TSS.percentile"
# [11] "gex_barcode"                       "atac_barcode"
# [13] "is_cell"                           "excluded_reason"
# [15] "gex_raw_reads"                     "gex_mapped_reads"
# [17] "gex_conf_intergenic_reads"         "gex_conf_exonic_reads"
# [19] "gex_conf_intronic_reads"           "gex_conf_exonic_unique_reads"
# [21] "gex_conf_exonic_antisense_reads"   "gex_conf_exonic_dup_reads"
# [23] "gex_exonic_umis"                   "gex_conf_intronic_unique_reads"
# [25] "gex_conf_intronic_antisense_reads" "gex_conf_intronic_dup_reads"
# [27] "gex_intronic_umis"                 "gex_conf_txomic_unique_reads"
# [29] "gex_umis_count"                    "gex_genes_count"
# [31] "atac_raw_reads"                    "atac_unmapped_reads"
# [33] "atac_lowmapq"                      "atac_dup_reads"
# [35] "atac_chimeric_reads"               "atac_mitochondrial_reads"
# [37] "atac_fragments"                    "atac_TSS_fragments"
# [39] "atac_peak_region_fragments"        "atac_peak_region_cutsites"
# [41] "pct_reads_in_peaks"                "blacklist_fraction"
# [43] "high.tss"                          "nucleosome_group"
# [45] "nCount_HTO"                        "nFeature_HTO"
# [47] "HTO_maxID"                         "HTO_secondID"
# [49] "HTO_margin"                        "HTO_classification"
# [51] "HTO_classification.global"         "hash.ID"
# [53] "nCount_mito"                       "nFeature_mito"
# [55] "mtDNA_depth"

VlnPlot(mtx, "mtDNA_depth", pt.size = 0.1) + scale_y_log10()
ggsave(filename="vlnplot_mtDNA_depth.jpg", path="/suffolk/Homes/ec474/scripts/mgatk/figures")

#We can look at the mitochondrial sequencing depth for each cell, and further subset the cells based on mitochondrial sequencing depth.

###################
##Calculate Heteroplasmy at mt.3243 allele position 
# as mutation is from A to G, sum all the G and the total, and make percentage

# suffolk[/home/ec474]% zcat mtx_mgatk_1.A.txt.gz | head -5
#position, CB, forward, reverse strand counts 
# 2,ATTGGTTCATGAAGTA-1,0,1
# 5,ATTGGTTCATGAAGTA-1,0,1
# 7,ATTGGTTCATGAAGTA-1,0,1
# 13,ATTGGTTCATGAAGTA-1,1,1
# 16,ATTGGTTCATGAAGTA-1,1,1



callHeteroplasmyByPosition <- function(seurat.obj, position, depth_cutoff) {
       # call heteroplasmy at a given position

       DefaultAssay(seurat.obj) <- "mito"

       mito_df <- as.data.frame(t(as.data.frame(seurat.obj@assays$mito@counts[grep(paste0("-", position, "-"), rownames(seurat.obj@assays$mito@counts)),])))
       mito_df$depth <- rowSums(mito_df)
       mito_df$allele_A <- mito_df[,1] + mito_df[,5] # sum forwards and reverse strands
       mito_df$allele_G <- mito_df[,4] + mito_df[,8]
       mito_df$Heteroplasmy <- mito_df$allele_G / (mito_df$allele_G + mito_df$allele_A) # MUT / (WT + MUT)
       mito_df$Heteroplasmy <- ifelse(mito_df$depth >= depth_cutoff, mito_df$Heteroplasmy, NA)

	return(mito_df)
}

# The function specifies -3243- so that grep will not pick up for instance numbers that contain this 
# but have additional numbers either before or after, such as 13243, or 32430

#mtx@assays$mito@counts[grep(3243, rownames(mtx@assays$mito@counts)),]
#16 x 221 sparse Matrix of class "dgCMatrix"
# This means that the matrix is subset by all those rows with position that includes "3243"
# Output
# A-3243-fwd  3 5 2 10 . 7 8 10 4 2 3 3  8 1 6  8 13 2 6 2 2 7 3 4 5 5 3 5 6 2 5
# A-13243-fwd . . .  . . . .  . . . . .  . . .  .  . . . . . . . . . . . . . . .
# C-3243-fwd  . . .  . . . .  . . . . .  . . .  .  . . . . . . . . . . . . . . .
# C-13243-fwd . . .  . . . .  . . . . .  . . .  .  . . . . . . . . . . . . . . .
# T-3243-fwd  . . .  . . . .  . . . . .  . . .  .  . . . . . . . . . . . . . . .
# T-13243-fwd . . .  . . . .  . . 0 . .  . . .  .  1 . . . . . . . . . . . . . .
# G-3243-fwd  . . .  . . . .  . . . . .  . . .  .  . . . . . . . . . . . . . . .
# G-13243-fwd 1 5 4  7 5 2 3  8 6 4 1 3  8 2 3 10  3 2 2 . 2 9 3 2 6 2 4 3 7 6 3

# head(colnames(mtx@assays$mito@counts))
# [1] "AAACGCGCAGCACGTT-1" "AAACGCGCAGGCAAGC-1" "AACAGGATCGCTAGAT-1"
# [4] "AACCGCTCAACGTGCT-1" "AACTACTCAGGAACCA-1" "AAGAACAGTGTTTCAC-1"

# mito.df <- as.data.frame(t(as.data.frame(mtx@assays$mito@counts[grep(-3243-, rownames(mtx@assays$mito@counts)),])))
#                       A-3243-fwd C-3243-fwd T-3243-fwd G-3243-fwd A-3243-rev
# AAACGCGCAGCACGTT-1          3          0          0          0          3
# AAACGCGCAGGCAAGC-1          5          0          0          0          3
# AACAGGATCGCTAGAT-1          2          0          0          0          2
# AACCGCTCAACGTGCT-1         10          0          0          0          9
# AACTACTCAGGAACCA-1          0          0          0          0          0
# AAGAACAGTGTTTCAC-1          7          0          0          0          4
#                    C-3243-rev T-3243-rev G-3243-rev
# AAACGCGCAGCACGTT-1          0          0          0
# AAACGCGCAGGCAAGC-1          0          0          0
# AACAGGATCGCTAGAT-1          0          0          0
# AACCGCTCAACGTGCT-1          0          0          0
# AACTACTCAGGAACCA-1          0          0          0
# AAGAACAGTGTTTCAC-1          0          0          0

#Therefore columns of interests are A-foward and A-reverse (complementary bases and strands) 
# And G-forward and G-reverse

#Therefore, rowSums(mito.df) would give the total count of all bases at that position for that particular cell
# length(mito.df$depth)
# [1] 221


het <- callHeteroplasmyByPosition(mtx, 3243)
colnames(het)
#  [1] "A-3243-fwd"   "C-3243-fwd"   "T-3243-fwd"   "G-3243-fwd"   "A-3243-rev"
#  [6] "C-3243-rev"   "T-3243-rev"   "G-3243-rev"   "depth"        "allele_A"
# [11] "allele_G"     "Heteroplasmy"

head(het) #heteroplasmy per cell
#   ...                Heteroplasmy 
# TCCAGGATCATCCTAT-1   0.00000000
# TCCAGGTCAATCTCTC-1   0.00000000
# TCTAAGGGTGCGCGTA-1   0.05555556
# TGAGCAAAGCTTATGA-1   0.00000000
# TGCACCTTCGATTTAG-1   0.00000000
# TGCAGGCTCTCCATGC-1   0.00000000


colSums(het)
  # A-3243-fwd   C-3243-fwd   T-3243-fwd   G-3243-fwd   A-3243-rev   C-3243-rev
  #       1178            0            0            9          712            0
  # T-3243-rev   G-3243-rev        depth     allele_A     allele_G Heteroplasmy
  #          0            1         1900         1890           10           NA



#Plot depth, and optimise the depth cut off, violn or density plot
#Published literature set this btw 10-20

mtx <- AddMetaData(mtx, metadata = het$Heteroplasmy, col.name = "Heteroplasmy")
#Check CB are in the same order

range(het$Heteroplasmy, na.rm=TRUE)
# [1] 0.0000000 0.4444444

library(ggplot2)
VlnPlot(mtx, "Heteroplasmy", pt.size = 0.1) + labs(title="scHeteroplasmy") + scale_y_continuous(labels = scales::percent)
ggsave(filename="vlnplot_scHeteroplasmy_depth_5.jpg", path="/suffolk/Homes/ec474/scripts/mgatk/figures")


#Average heteroplasmy of sample
het_sum <- mean(het$Heteroplasmy, na.rm=TRUE)
cells <- length(het$Heteroplasmy)
het_av <- het_sum / cells


#mtx <- saveRDS(file= "/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000_singlet_HASH1_RNA_ATAC_QC.rds")






#DO NOT SUBSET BY MEAN MITOCHONDRIAL DEPTH AS YOU ONLY CARE ABOUT THE DEPTH AT THE PARTICULAR MUTATED SITE
# filter cells based on mitochondrial depth
mtx <- subset(mtx, mtDNA_depth >= 10) #Only 47 samples left out of 221
mtx
# An object of class Seurat
# 267363 features across 47 samples within 4 assays
# Active assay: ATAC (98206 features, 0 variable features)
#  2 layers present: counts, data
#  3 other assays present: RNA, HTO, mito

#Dimensionality reduction and clustering of ATAC assay 
mtx <- RunTFIDF(mtx)
# Warning message:
# In RunTFIDF.default(object = GetAssayData(object = object, slot = "counts"),  :
#   Some features contain 0 total counts

mtx <- FindTopFeatures(mtx, min.cutoff = 10)
mtx <- RunSVD(mtx)
# Running SVD
# Warning in irlba(A = t(x = object), nv = n, work = irlba.work, tol = tol) :
#   You're computing too large a percentage of total singular values, use a standard svd instead.
# Scaling cell embeddings

mtx <- RunUMAP(mtx, reduction = "lsi", dims = 2:50)
# Error in Embeddings(object[[reduction]])[, dims] :
#   subscript out of bounds

# Selected dims = 2:30 instead
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# 12:26:52 UMAP embedding parameters a = 0.9922 b = 1.112
# 12:26:52 Read 47 rows and found 29 numeric columns
# 12:26:52 Using Annoy for neighbor search, n_neighbors = 30
# 12:26:52 Building Annoy index with metric = cosine, n_trees = 50
# 12:26:52 Writing NN index file to temp file /tmp/Rtmp2SuoXu/file3dff3d8d5344
# 12:26:52 Searching Annoy index using 1 thread, search_k = 3000
# 12:26:52 Annoy recall = 100%
# 12:26:54 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 12:26:56 Initializing from normalized Laplacian + noise (using RSpectra)
# 12:26:56 Commencing optimization for 500 epochs, with 1428 positive edges
# Using method 'umap'
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 12:26:57 Optimization finished

mtx <- FindNeighbors(mtx, reduction = "lsi", dims = 2:30)
# Computing nearest neighbor graph
# Computing SNN

mtx <- FindClusters(mtx, resolution = 0.7, algorithm = 3)
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# Number of nodes: 47
# Number of edges: 1081

# Running smart local moving algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Maximum modularity in 10 random starts: 0.3000
# Number of communities: 1
# Elapsed time: 0 seconds

DimPlot(mtx, label = TRUE) + NoLegend()
ggsave(filename="umap_atac.jpg", path="/suffolk/Homes/ec474/scripts/mgatk/figures")
#There are no cell clusters in this UMAP


#Generate gene scores 
#To help interpret these clusters of cells, and assign a cell type label, we’ll estimate gene activities by summing the DNA accessibility in the gene body and promoter region.
# compute gene accessibility
gene.activities <- GeneActivity(mtx)

# add to the Seurat object as a new assay
mtx[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)

# An object of class Seurat
# 286970 features across 47 samples within 5 assays
# Active assay: ATAC (98206 features, 12215 variable features)
#  2 layers present: counts, data
#  4 other assays present: RNA, HTO, mito, GeneActivity
#  2 dimensional reductions calculated: lsi, umap

mtx <- NormalizeData(
  object = mtx,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(mtx$nCount_RNA)
)

#Visualise intere
    # EPCAM is a marker for epithelial cells
    # TREM1 is a meyloid marker
    # PTPRC = CD45 is a pan-immune cell marker
    # IL1RL1 is a basophil marker
    # GATA3 is a Tcell maker

DefaultAssay(mtx) <- 'GeneActivity'

FeaturePlot(
  object = mtx,
  features = c('TREM1', 'EPCAM', "PTPRC", "IL1RL1","GATA3", "KIT"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 2
)
ggsave(filename="gene_activity_plot.jpg", path="/suffolk/Homes/ec474/scripts/mgatk/figures")


# I do not have clusters --- so no point doing the steps below
# Using these gene score values, we can assign cluster identities:
mtx <- RenameIdents(
  object = mtx,
  '0' = 'Epithelial',
  '1' = 'Epithelial',
  '2' = 'Basophil',
  '3' = 'Myeloid_1',
  '4' = 'Myeloid_2',
  '5' = 'Tcell'
)

# One of the myeloid clusters has a lower percentage of fragments in peaks, 
# as well as a lower overall mitochondrial sequencing depth and a different nucleosome banding pattern.

p1 <- FeatureScatter(crc, "mtDNA_depth", "pct_reads_in_peaks") + ggtitle("") + scale_x_log10()
p2 <- FeatureScatter(crc, "mtDNA_depth", "nucleosome_signal") + ggtitle("") + scale_x_log10()

p1 + p2 + plot_layout(guides = 'collect')

#Find informative mtDNA variants
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = mito.data$refallele)
VariantPlot(variants = variable.sites)

# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

high.conf[,c(1,2,5)]

# Compute the variant allele frequency for each cell
crc <- AlleleFreq(
  object = crc,
  variants = high.conf$variant,
  assay = "mito"
)
crc[["alleles"]]

# Visualize the variants
DefaultAssay(crc) <- "alleles"
alleles.view <- c("12889G>A", "16147C>T", "9728C>T", "9804G>A")
FeaturePlot(
  object = crc,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 4
) & NoLegend()

DoHeatmap(crc, features = rownames(crc), slot = "data", disp.max = 1) +
  scale_fill_viridis_c()

# Identifying clones
DefaultAssay(tf1) <- "alleles"
tf1 <- FindClonotypes(tf1)

table(Idents(tf1))

DoHeatmap(tf1, features = VariableFeatures(tf1), slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c()

# Find differentially accessible peaks between clones
DefaultAssay(tf1) <- "peaks"

# find peaks specific to one clone
markers.fast <- FoldChange(tf1, ident.1 = 2)
markers.fast <- markers.fast[order(markers.fast$avg_log2FC, decreasing = TRUE), ] # sort by fold change
head(markers.fast)

CoveragePlot(
  object = tf1,
  region = rownames(markers.fast)[1],
  extend.upstream = 2000,
  extend.downstream = 2000
)





