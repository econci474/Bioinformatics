#Running the QC on the ATAC slot for the Seurat object that already had RNA QC done to it

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(hdf5r)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

###############
mtx <- readRDS("/suffolk/RawGenomicsF/ec474/mitox_sce/mtx_qc.rds")

mtx
#An object of class Seurat
#134807 features across 9907 samples within 2 assays
#Active assay: RNA (36601 features, 0 variable features)
# 1 layer present: counts
# 1 other assay present: ATAC

DefaultAssay(mtx) <- "ATAC"
mtx
#An object of class Seurat
#134807 features across 9907 samples within 2 assays
#Active assay: ATAC (98206 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA

metadata <- read.csv(file="/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/per_barcode_metrics.csv", header = TRUE,
  row.names = 1)

head(colnames(metadata))
#[1] "gex_barcode"      "atac_barcode"     "is_cell"          "excluded_reason"
#[5] "gex_raw_reads"    "gex_mapped_reads"

mtx <- AddMetaData(mtx, metadata)

mtx[['ATAC']]
#ChromatinAssay data with 98206 features for 9907 cells
#Variable features: 0
#Genome: hg38
#Annotation present: TRUE
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
#Found 9907 cell barcodes
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
colnames(mtx[[]])
#[1] "orig.ident"                        "nCount_RNA"
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

#variables in cell-ranger-arc are different, see 
#https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics
 
mtx$pct_reads_in_peaks <- mtx$atac_peak_region_fragments / mtx$atac_fragments * 100

mtx$blacklist_fraction <- FractionCountsInRegion(
  object = mtx, 
  assay = 'ATAC',
  regions = blacklist_hg38_unified
)

#Visualise variables stored in object metadata
#Set quantiles to find suitable cutoff values for QC 

DensityScatter(mtx, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(filename="Density_nCount_ATAC_TSS.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

#Inspect TSS enrichment scores by grouping cells based on score and plotting accessibility signal over all TSS sites

mtx$high.tss <- ifelse(mtx$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(mtx, group.by = 'high.tss') + NoLegend()
ggsave(filename="TSSplot.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

#fragment length periodicity for all the cells, and group by cells with high or low nucleosomal signal strength. 
#cells that are outliers for the mononucleosomal / nucleosome-free ratio (based on the plots above) have different nucleosomal banding patterns. 
#The remaining cells exhibit a pattern that is typical for a successful ATAC-seq experiment.

#Fragment size distribution typically reflects nucleosome binding pattern showing enrichment around values corresponding to fragments bound to a single nucleosome (between 147 bp and 294 bp) 
#as well as nucleosome-free fragments (shorter than 147 bp). 
#The ratio of mono-nucleosome cut fragments to nucleosome-free fragments can be called nucleosome signal, and it can be estimated using a subset of fragments with 

frag <- FragmentHistogram(object = mtx,assay='ATAC') + geom_vline(xintercept=147, color="red")
frag + scale_fill_manual(values="lightblue")
ggsave(filename="Fragment_histogram_total.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

mtx$nucleosome_group <- ifelse(mtx$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = mtx, group.by = 'nucleosome_group')
ggsave(filename="Fragment_histogram_nucleosome_group.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

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
#  9907

library(patchwork)
vln <- VlnPlot(
  object = mtx,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.01,
  ncol = 5
)
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "Before Filtering")
ggsave(filename="Vln_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/", width=40,units="cm")

#without dots and with hline now
#nCount_ATAC, trip top and bottom 5% as per Density plot
VlnPlot(mtx, assay="ATAC", features="nCount_ATAC", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=c(100000,1632), color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nCount_ATAC_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

sum(mtx@meta.data$nCount_ATAC>100000)
#[1] 74
sum(mtx@meta.data$nCount_ATAC>52849, mtx@meta.data$nCount_ATAC<1632)
#[1] 992
sum(mtx@meta.data$nCount_ATAC<1632) #bottom 5%
#[1] 496

#TSS.enrichment
VlnPlot(mtx, assay="ATAC", features="TSS.enrichment", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=2.5, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_TSS_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

sum(mtx@meta.data$TSS.enrichment<3)
#[1] 56

sum(mtx@meta.data$TSS.enrichment<2.5)
#[1] 11

#blacklist_fraction
VlnPlot(mtx, assay="ATAC", features="blacklist_fraction", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=0.005, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_blacklist_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

sum(mtx@meta.data$blacklist_fraction>0.005)
#[1] 27

#nucleosome_signal
VlnPlot(mtx, assay="ATAC", features="nucleosome_signal", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=2, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nucleosome_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

sum(mtx@meta.data$nucleosome_signal>2)
#[1] 2

#pct_reads_in_peaks
VlnPlot(mtx, assay="ATAC", features="pct_reads_in_peaks", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=20, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_pct_reads_in_peaks_before.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

sum(mtx@meta.data$pct_reads_in_peaks<20)
#[1] 77

#filtering 
mtx <- subset(
  x = mtx,
  subset = nCount_ATAC > 1632 &
    nCount_ATAC < 100000 &
    pct_reads_in_peaks > 20 &
    blacklist_fraction < 0.005 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2.5
)
mtx
#An object of class Seurat
#98239 features across 9060 samples within 1 assay
#Active assay: peaks (98239 features, 0 variable features)
# 2 layers present: counts, data

vln <- VlnPlot(
  object = mtx,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.01,
  ncol = 5
)
vln <- vln & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
vln + NoLegend() + plot_annotation(title = "After Filtering")
ggsave(filename="Vln_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/", width=40,units="cm")

#after filtering violins
#nCount_peaks
VlnPlot(mtx, assay="ATAC", features="nCount_ATAC", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=c(100000,1632), color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nCount_ATAC_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

#TSS.enrichment
VlnPlot(mtx, assay="ATAC", features="TSS.enrichment", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=2.5, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_TSS_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

#blacklist_fraction
VlnPlot(mtx, assay="ATAC", features="blacklist_fraction", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=0.005, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_blacklist_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

#nucleosome_signal
VlnPlot(mtx, assay="ATAC", features="nucleosome_signal", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=4, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_nucleosome_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

#pct_reads_in_peaks
VlnPlot(mtx, assay="ATAC", features="pct_reads_in_peaks", pt.size = 0.01) + NoLegend() + geom_hline(yintercept=20, color="red") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("VlnPlot_pct_reads_in_peaks_after.jpg", path="/suffolk/Homes/ec474/scripts/seurat_RNA_ATAC_QC/figures_atac_integrated_rna/")

saveRDS(mtx, file="/suffolk/WorkGenomicsD/ec474/mitox_qc/mtx_rna_atac_qc.rds")
