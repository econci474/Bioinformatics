#!/bin/bash

#SBATCH -D /suffolk/WorkGenomicsD/ec474/mgatk/
#SBATCH -J MGATK
#SBATCH -o run_mgatk.log           #Job output
#SBATCH -t 300:00:00         #Max wall time for entire job
#SBATCH --partition=sledges7  #this could be any other partition on the cluster
#SBATCH --mem=50M #with this you specify the memory per node
#SBATCH -N 1 #number of nodes. Note all your jobs could use a max of 4 nodes on our MBU cluster.
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user ec474@mrc-mbu.cam.ac.uk  #your email

# Activate the virtual environment
source /home/ec474/venv3/bin/activate
pip list

input_dir="/suffolk/WorkGenomicsD/ec474/mitox_count_masked"
output_dir="/suffolk/WorkGenomicsD/ec474/mgatk"

mgatk tenx -i ${input_dir}/outs/atac_possorted_bam.bam \
-n mtx_mgatk_1 -o ${output_dir} \
-c 12 \
-bt CB \
-b ${input_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv

echo DONE!