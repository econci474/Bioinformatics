#Following Signac vignette: https://stuartlab.org/signac/articles/pbmc_vignette
#Running QC on ATAC slot only without any RNA processing to learn how to do it

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(hdf5r)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

#mtx <- readRDS("/suffolk/RawGenomicsF/ec474/mitox_sce/mtx.rds")

#Try from beginning rather than loading in object, in case of any errors
counts <- Read10X_h5("/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/filtered_feature_bc_matrix.h5", use.names=TRUE, unique.features = TRUE)

mtx.rna <- counts$`Gene Expression`
mtx.atac <- counts$Peaks

mtx <- CreateSeuratObject(counts = mtx.rna, assay="RNA")


###############

#fragment file
fragpath = "/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/atac_fragments.tsv.gz"

#for the multiome kit the metadata can be found in the cellranger-arc output: outs/per_barcode_metrics.csv

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = fragpath,
  min.cells = 10,
  min.features = 200
)
#Error in `separate()`:
#! Can't extract columns that don't exist.
#✖ Column `ranges` doesn't exist.
#Run `rlang::last_trace()` to see where the error occurred.

chrom_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  min.cells = 10,
  min.features = 200
)
#Computing hash
#Checking for 10055 cell barcodes

metadata <- read.csv(file="/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/per_barcode_metrics.csv", header = TRUE,
  row.names = 1)

head(colnames(metadata))
#[1] "gex_barcode"      "atac_barcode"     "is_cell"          "excluded_reason"
#[5] "gex_raw_reads"    "gex_mapped_reads"

mtx <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

mtx
#An object of class Seurat
#98239 features across 10055 samples within 1 assay
#Active assay: peaks (98239 features, 0 variable features)
# 2 layers present: counts, data

mtx[['peaks']]
#ChromatinAssay data with 98239 features for 10055 cells
#Variable features: 0
#Genome:
#Annotation present: FALSE
#Motifs present: FALSE
#Fragment files: 1

granges(mtx)
#GRanges object with 98239 ranges and 0 metadata columns:
#            seqnames        ranges strand
#               <Rle>     <IRanges>  <Rle>
#      [1]       chr1    9803-10701      *
#      [2]       chr1 181010-181757      *
#      [3]       chr1 778288-779195      *
#      [4]       chr1 816891-817776      *
#      [5]       chr1 819983-820697      *
#      ...        ...           ...    ...
#  [98235] KI270713.1     2759-3531      *
#  [98236] KI270713.1     3949-4851      *
#  [98237] KI270713.1   21459-22374      *
#  [98238] KI270713.1   26041-26910      *
#  [98239] KI270713.1   36931-37819      *
#  -------
#  seqinfo: 33 sequences from an unspecified genome; no seqlengths

# extract gene annotations from EnsDb
annotation <- GetGRangesFromEnsDb(ensdb= EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
seqlevels(annotation) <- paste0('chr', seqlevels(annotation)) 
genome(annotation) <- "hg38"
Annotation(mtx) <- annotation

Annotation(mtx)
#GRanges object with 3021151 ranges and 5 metadata columns:
#                  seqnames        ranges strand |           tx_id   gene_name
#                     <Rle>     <IRanges>  <Rle> |     <character> <character>
#  ENSE00001489430     chrX 276322-276394      + | ENST00000399012      PLCXD1
#  ENSE00001536003     chrX 276324-276394      + | ENST00000484611      PLCXD1
#  ENSE00002160563     chrX 276353-276394      + | ENST00000430923      PLCXD1
#  ENSE00001750899     chrX 281055-281121      + | ENST00000445062      PLCXD1
#  ENSE00001489388     chrX 281192-281684      + | ENST00000381657      PLCXD1
#              ...      ...           ...    ... .             ...         ...
#  ENST00000361739    chrMT     7586-8269      + | ENST00000361739      MT-CO2
#  ENST00000361789    chrMT   14747-15887      + | ENST00000361789      MT-CYB
#  ENST00000361851    chrMT     8366-8572      + | ENST00000361851     MT-ATP8
#  ENST00000361899    chrMT     8527-9207      + | ENST00000361899     MT-ATP6
#  ENST00000362079    chrMT     9207-9990      + | ENST00000362079      MT-CO3
#                          gene_id   gene_biotype     type
#                      <character>    <character> <factor>
#  ENSE00001489430 ENSG00000182378 protein_coding     exon
#  ENSE00001536003 ENSG00000182378 protein_coding     exon
#  ENSE00002160563 ENSG00000182378 protein_coding     exon
#  ENSE00001750899 ENSG00000182378 protein_coding     exon
#  ENSE00001489388 ENSG00000182378 protein_coding     exon
#              ...             ...            ...      ...
#  ENST00000361739 ENSG00000198712 protein_coding      cds
#  ENST00000361789 ENSG00000198727 protein_coding      cds
#  ENST00000361851 ENSG00000228253 protein_coding      cds
#  ENST00000361899 ENSG00000198899 protein_coding      cds
#  ENST00000362079 ENSG00000198938 protein_coding      cds
#  -------
#  seqinfo: 25 sequences (1 circular) from hg38 genome


#Nucleosome banding pattern: 
#The histogram of DNA fragment sizes (determined from the paired-end sequencing reads) should exhibit a strong nucleosome banding pattern corresponding to the length of DNA wrapped around a single nucleosome. 
#We calculate this per single cell, and quantify the approximate ratio of mononucleosomal to nucleosome-free fragments (stored as nucleosome_signal)

mtx <- NucleosomeSignal(object = mtx)
#Found 10055 cell barcodes
#Done Processing 50 million lines

#Transcriptional start site (TSS) enrichment score. 
#The ENCODE project has defined an ATAC-seq targeting score based on the ratio of fragments centered at the TSS to fragments in TSS-flanking regions 
#(see https://www.encodeproject.org/data-standards/terms/). 
#Poor ATAC-seq experiments typically will have a low TSS enrichment score. 
#We can compute this metric for each cell with the TSSEnrichment() function, and the results are stored in metadata under the column name TSS.enrichment.

#important that fast=FALSE if you want to plot TSSplot
mtx <- TSSEnrichment(object = mtx, fast = FALSE)
#Extracting TSS positions
#Finding + strand cut sites
#Finding - strand cut sites
#Computing mean insertion frequency in flanking regions
#Normalizing TSS score

#Total number of fragments in peaks: A measure of cellular sequencing depth / complexity. 
#Cells with very few reads may need to be excluded due to low sequencing depth. 
#Cells with extremely high levels may represent doublets, nuclei clumps, or other artefacts.

#Fraction of fragments in peaks: 
#Represents the fraction of all fragments that fall within ATAC-seq peaks. 
#Cells with low values (i.e. <15-20%) often represent low-quality cells or technical artifacts that should be removed. 
#Note that this value can be sensitive to the set of peaks used.

#Ratio reads in genomic blacklist regions 
#The ENCODE project has provided a list of blacklist regions, representing reads which are often associated with artefactual signal. 
#Cells with a high proportion of reads mapping to these areas (compared to reads mapping to peaks) often represent technical artifacts and should be removed. 
#ENCODE blacklist regions for human (hg19 and GRCh38), mouse (mm10), Drosophila (dm3), and C. elegans (ce10) are included in the Signac package.

str(mtx)
head(mtx[[]])

#variables in cell-ranger-arc are different, see 
#https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics

mtx$pct_reads_in_peaks <- mtx$atac_peak_region_fragments / mtx$atac_fragments * 100

mtx$blacklist_fraction <- FractionCountsInRegion(
  object = mtx, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

#Visualise variables stored in object metadata
#Set quantiles to find suitable cutoff values for QC 

DensityScatter(mtx, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(filename="Density_nCount_peaks_TSS.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#Inspect TSS enrichment scores by grouping cells based on score and plotting accessibility signal over all TSS sites

mtx$high.tss <- ifelse(mtx$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(mtx, group.by = 'high.tss') + NoLegend()
ggsave(filename="TSSplot.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength. 
#cells that are outliers for the mononucleosomal / nucleosome-free ratio (based on the plots above) have different nucleosomal banding patterns. 
#The remaining cells exhibit a pattern that is typical for a successful ATAC-seq experiment.

mtx$nucleosome_group <- ifelse(mtx$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = mtx, group.by = 'nucleosome_group')
ggsave(filename="Fragment_histogram.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#Warning messages:
#1: Removed 5066 rows containing non-finite outside the scale range (`stat_bin()`).
#2: Removed 2 rows containing missing values or values outside the scale range (`geom_bar()`).
#Saving 7 x 7 in image
#Warning messages:
#1: Removed 5066 rows containing non-finite outside the scale range (`stat_bin()`).
#2: Removed 2 rows containing missing values or values outside the scale range (`geom_bar()`).

#Scores represent hybridization strength and are inversely proportional to the amount of cleavage; 
#i.e., high-scoring fragments are nucleosome forming and low-scoring fragments are non-nucleosome forming.

#This Looks problematic!
table(mtx$nucleosome_group)
#NS < 4
# 10055

library(patchwork)
vln <- VlnPlot(
  object = mtx,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "Before Filtering")
ggsave(filename="Vln_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/", width=40,units="cm")

#without dots and with hline now
#nCount_peaks
VlnPlot(mtx, assay="peaks", features="nCount_peaks", pt.size = 0) + NoLegend() + geom_hline(yintercept=c(100000,3000), color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nCount_peaks_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#TSS.enrichment
VlnPlot(mtx, assay="peaks", features="TSS.enrichment", pt.size = 0) + NoLegend() + geom_hline(yintercept=3, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_TSS_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#blacklist_fraction
VlnPlot(mtx, assay="peaks", features="blacklist_fraction", pt.size = 0) + NoLegend() + geom_hline(yintercept=0.005, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_blacklist_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#nucleosome_signal
VlnPlot(mtx, assay="peaks", features="nucleosome_signal", pt.size = 0) + NoLegend() + geom_hline(yintercept=4, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nucleosome_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#pct_reads_in_peaks
VlnPlot(mtx, assay="peaks", features="pct_reads_in_peaks", pt.size = 0) + NoLegend() + geom_hline(yintercept=20, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_pct_reads_in_peaks_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#filtering 
mtx <- subset(
  x = mtx,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 20 &
    blacklist_fraction < 0.005 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
mtx
#An object of class Seurat
#98239 features across 9060 samples within 1 assay
#Active assay: peaks (98239 features, 0 variable features)
# 2 layers present: counts, data

vln <- VlnPlot(
  object = mtx,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "After Filtering")
ggsave(filename="Vln_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/", width=40,units="cm")

#after filtering violins
#nCount_peaks
VlnPlot(mtx, assay="peaks", features="nCount_peaks", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=c(100000,3000), color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nCount_peaks_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#TSS.enrichment
VlnPlot(mtx, assay="peaks", features="TSS.enrichment", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=3, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_TSS_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#blacklist_fraction
VlnPlot(mtx, assay="peaks", features="blacklist_fraction", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=0.005, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_blacklist_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#nucleosome_signal
VlnPlot(mtx, assay="peaks", features="nucleosome_signal", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=4, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nucleosome_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

#pct_reads_in_peaks
VlnPlot(mtx, assay="peaks", features="pct_reads_in_peaks", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=20, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_pct_reads_in_peaks_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac/")

saveRDS(mtx, file="/suffolk/WorkGenomicsD/ec474/mitox_qc/mtx_atac_qc.rds")


#Normalisation and linear dimensionality reduction 
mtx <- RunTFIDF(mtx)
mtx <- FindTopFeatures(mtx, min.cutoff = 'q0')
mtx <- RunSVD(mtx)

#The first LSI component often captures sequencing depth (technical variation) rather than biological variation. If this is the case, the component should be removed from downstream analysis. 
#We can assess the correlation between each LSI component and sequencing depth using the DepthCor() function:
DepthCor(mtx)
ggsave(filename="depthcorrelation.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_dimred/")

head(mtx@reductions)

#Can't elbow plot
ElbowPlot(mtx, ndims=50, reduction="lsi")
ggsave(filename="elbowplot.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_dimred/")


#
mtx <- RunUMAP(object = mtx, reduction = 'lsi', dims = 2:30)
mtx <- FindNeighbors(object = mtx, reduction = 'lsi', dims = 2:30)
mtx <- FindClusters(object = mtx, verbose = FALSE, algorithm = 3)
DimPlot(object = mtx, label = TRUE) + NoLegend() + labs(title="UMAP of ATAC LSI")
ggsave(filename="umap_atac.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_dimred/")

#annotating and interpreting clusters is more challenging in scATAC-seq data as much less is known about the functional roles of noncoding genomic regions than is known about protein coding regions (genes).
#We can try to quantify the activity of each gene in the genome by assessing the chromatin accessibility associated with the gene, and create a new gene activity assay derived from the scATAC-seq data. 
#summing the fragments intersecting the gene body and promoter region 
#(we also recommend exploring the Cicero tool, which can accomplish a similar goal, and we provide a vignette showing how to run Cicero within a Signac workflow here).

#To create a gene activity matrix, we extract gene coordinates and extend them to include the 2 kb upstream region (as promoter accessibility is often correlated with gene expression). 
#We then count the number of fragments for each cell that map to each of these regions, using the using the FeatureMatrix() function. 
#These steps are automatically performed by the GeneActivity() function:

#Gene activity matrix
gene.activities <- GeneActivity(mtx)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
mtx[['RNA']] <- CreateAssayObject(counts = gene.activities)
mtx <- NormalizeData(
  object = mtx,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(mtx$nCount_RNA)
)

DefaultAssay(mtx) <- 'RNA'
FeaturePlot(
  object = mtx,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
ggsave(filename="featureplot_geneactivities_atac.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_dimred/", width=30, units="cm")

saveRDS(mtx, file="/suffolk/WorkGenomicsD/ec474/mitox_qc/mtx_atac_qc_dimred_multiplexed.rds")

#metadata to compare across with demultiplexed object
DefaultAssay(mtx) <- "peaks"
#An object of class Seurat
#117846 features across 9060 samples within 2 assays
#Active assay: peaks (98239 features, 98239 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: lsi, umap

head(colnames(mtx))
#[1] "AAACAGCCACCTAATG-1" "AAACAGCCATAGACCC-1" "AAACAGCCATTAAAGG-1"
#[4] "AAACATGCAAGGTGGC-1" "AAACATGCAATCCTGA-1" "AAACATGCACCTCACC-1"

head(rownames(mtx))
#[1] "chr1-9803-10701"    "chr1-181010-181757" "chr1-778288-779195"
#[4] "chr1-816891-817776" "chr1-819983-820697" "chr1-822687-823558"

head(mtx@meta.data)
#[1] "orig.ident"                        "nCount_peaks"
# [3] "nFeature_peaks"                    "gex_barcode"
# [5] "atac_barcode"                      "is_cell"
# [7] "excluded_reason"                   "gex_raw_reads"
# [9] "gex_mapped_reads"                  "gex_conf_intergenic_reads"
#[11] "gex_conf_exonic_reads"             "gex_conf_intronic_reads"
#[13] "gex_conf_exonic_unique_reads"      "gex_conf_exonic_antisense_reads"
#[15] "gex_conf_exonic_dup_reads"         "gex_exonic_umis"
#[17] "gex_conf_intronic_unique_reads"    "gex_conf_intronic_antisense_reads"
#[19] "gex_conf_intronic_dup_reads"       "gex_intronic_umis"
#[21] "gex_conf_txomic_unique_reads"      "gex_umis_count"
#[23] "gex_genes_count"                   "atac_raw_reads"
#[25] "atac_unmapped_reads"               "atac_lowmapq"
#[27] "atac_dup_reads"                    "atac_chimeric_reads"
#[29] "atac_mitochondrial_reads"          "atac_fragments"
#[31] "atac_TSS_fragments"                "atac_peak_region_fragments"
#[33] "atac_peak_region_cutsites"         "nucleosome_signal"
#[35] "nucleosome_percentile"             "TSS.enrichment"
#[37] "TSS.percentile"                    "high.tss"
#[39] "nucleosome_group"                  "pct_reads_in_peaks"
#[41] "blacklist_fraction"                "peaks_snn_res.0.8"
#[43] "seurat_clusters"                   "nCount_RNA"
#[45] "nFeature_RNA"

head(mtx@meta.data$atac_barcode)
#[1] "ACAGCGGGTTACTAGC-1" "ACAGCGGGTAAAGGTT-1" "ACAGCGGGTCATTACA-1"
#[4] "CATTTAGGTGCTAGCC-1" "CATTTAGGTAGACAAT-1" "CATTTAGGTATTGTCA-1"