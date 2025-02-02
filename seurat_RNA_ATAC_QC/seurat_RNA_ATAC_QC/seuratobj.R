
library(dplyr)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(hdf5r)
library(GenomicRanges)

#Create Seurat object containing RNA data
mtx.inputs <- Read10X_h5("/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/filtered_feature_bc_matrix.h5", use.names=TRUE, unique.features = TRUE)
mtx.rna <- mtx.inputs$`Gene Expression`
mtx.atac <- mtx.inputs$Peaks

mtx <- CreateSeuratObject(counts = mtx.rna, assay="RNA")

#fragment file
fragpath = "/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/atac_fragments.tsv.gz"

#only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(mtx.atac), sep = c(":", "-"))
seqlevels(grange.counts, pruning.mode = "coarse") <- standardChromosomes(grange.counts)
grange.use <- seqlevels(grange.counts)
chromosomes <- sapply(strsplit(rownames(mtx.atac), split=":"), `[`, 1)
mtx.atac <- mtx.atac[chromosomes %in% grange.use, ]

#Get annotations for hg38 with NUMT masking
annotation <- GetGRangesFromEnsDb(ensdb= EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation)) 
genome(annotation) <- "hg38"

Annotation(mtx.atac) <- annotation

chrom_assay <- CreateChromatinAssay(
    counts=mtx.atac, 
    sep = c(":", "-"), 
    fragments = fragpath, 
    annotation = annotation)

genome(chrom_assay) <- "hg38"

#add metadata to seurat object from cell ranger arc: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics 
metadata <- read.csv(file="/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/per_barcode_metrics.csv", header = TRUE,
  row.names = 1)

#Is this vadlid??
LayerData(chrom_assay, layer = "metadata") <- metadata

mtx.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = metadata
)

#Is this valid?? 
mtx[["ATAC"]] <- mtx.atac[["ATAC"]]

#Add seurat object to another seurat object
mtx[["ATAC"]] <- chrom_assay

#or Add Meta data to seurat object later? 
mtx <- AddMetaData(mtx, metadata)

Assays(mtx)
saveRDS(mtx, file = "/suffolk/RawGenomicsF/ec474/mitox_sce/mtx.rds")
