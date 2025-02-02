
library(dplyr)
library(ggplot2)
library(Seurat)

#New object, try to subset the metadata only (this should still be accurate?) using dplyr package
mtx_umis_clara_0.95_1000 <- readRDS("/suffolk/WorkGenomicsD/ec474/demultiplex/mtx_umis_clara_0.95_1000.rds")


#Subset HASH-1 patient
df <- as.data.frame(mtx_umis_clara_0.95_1000@meta.data)
df_subset <- df %>% dplyr::filter(HTO_classification.global=="Singlet" & HTO_maxID=="NPC-HASH1")

df_subset_cell_id <- rownames(df_subset)
head(df_subset_cell_id)
# [1] "AAACGCGCAGCACGTT-1" "AAACGCGCAGGCAAGC-1" "AACAGGATCGCTAGAT-1"
# [4] "AACCGCTCAACGTGCT-1" "AACTACTCAGGAACCA-1" "AAGAACAGTGTTTCAC-1"
length(df_subset_cell_id)
# [1] 221

write.table(df_subset_cell_id, file="/suffolk/WorkGenomicsD/ec474/mgatk/hash_barcodes/hash1_barcodes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Subset HASH-2 patient
df <- mtx_umis_clara_0.95_1000@meta.data
df_subset <- df %>% dplyr::filter(HTO_classification.global=="Singlet" & HTO_maxID=="NPC-HASH2")

df_subset_cell_id <- rownames(df_subset)
head(df_subset_cell_id)
# [1] "AACAGATAGTTATCTC-1" "AACGGTAAGGTCGATT-1" "AACTAGTGTTGCCTCA-1"
# [4] "AAGCCTGTCCGTCCAT-1" "AAGTTAGCAGTATGTT-1" "AAGTTTGTCAGCAAAG-1"
length(df_subset_cell_id)
# [1] 159

write.table(df_subset_cell_id, file="/suffolk/WorkGenomicsD/ec474/mgatk/hash_barcodes/hash2_barcodes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
