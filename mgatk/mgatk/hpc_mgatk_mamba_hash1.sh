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
module load miniconda/2024-02-20
conda activate /home/ec474/miniforge3/envs/MGATK

# Print environment details
echo "Conda environment:"
conda info

echo "Installed packages:"
conda list

echo "Python version:"
python --version

echo "Python path:"
which python

echo "Pysam version:"
pysam --version

echo "Pysam path:"
which pysam

# Source JAVA package
JAVA_HOME=/export/MBU/software/anaconda3/pkgs/openjdk-8.0.121-1 
export JAVA_HOME
export PATH=$JAVA_HOME/bin:$PATH

#Code
input_dir=/suffolk/WorkGenomicsD/ec474/mitox_count_masked
output_dir=/suffolk/WorkGenomicsD/ec474/mgatk
barcodes=/suffolk/WorkGenomicsD/ec474/mgatk/hash_barcodes

mgatk tenx -i ${input_dir}/outs/atac_possorted_bam.bam \
-n mtx_mgatk_1 -o ${output_dir} \
-c 12 \
-bt CB \
-b ${barcodes}/hash1_barcodes.tsv

echo DONE!