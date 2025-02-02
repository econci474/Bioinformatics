#!/bin/bash

#SBATCH -D /suffolk/RawGenomicsF/ec474/cite-seq_count5/
#SBATCH -J CITE_SC
#SBATCH -o run_citeseqcount.log           #Job output
#SBATCH -t 300:00:00         #Max wall time for entire job
#SBATCH --partition=sledges7  #this could be any other partition on the cluster
#SBATCH --mem=50M #with this you specify the memory per node
#SBATCH -N 1 #number of nodes. Note all your jobs could use a max of 4 nodes on our MBU cluster.
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user ec474@mrc-mbu.cam.ac.uk  #your email

module load miniconda/2024-02-20
conda activate /home/ec474/miniforge3/envs/CITESEQCOUNT
conda list

bash /home/ec474/scripts/citeseqcount.sh

echo DONE!