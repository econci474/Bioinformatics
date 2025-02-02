input_dir="/suffolk/WorkGenomicsD/ec474/mitox_count_masked"
output_dir="/suffolk/WorkGenomicsD/ec474/mgatk"

mgatk tenx -i ${input_dir}/outs/atac_possorted_bam.bam \
-n mtx_mgatk_1 -o ${output_dir} \
-c 12 \
-bt CB \
-b ${input_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv

# -i is the filepath to input .bam file, which can be found in cellranger-arc ATAC outs folder
# -n output file name 
# -o output file directory
# -c is the number of cores used for genotyping
# - bt CB indicates that the CB SAM tag denotes the barcodes per single-cell deualt for 10X .bam files
# -ub UB  is added on if you have a scRNA seq bam file, which will require UMI-aware PCR deduplication. 
# Therefore we need to inform mgatk to look for error-corrected UMIs 
# -b are the barcodes a priori from the Cellranger knee-call 

