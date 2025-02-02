#!/bin/bash

#SBATCH -D /work/WorkGenomicsD/ec474
#SBATCH -J MTX
#SBATCH -o run_MTX_output.log           #Job output
#SBATCH -t 300:00:00         #Max wall time for entire job
#SBATCH --partition=sledges7  #this could be any other partition on the cluster
#SBATCH --mem=50M #with this you specify the memory per node
#SBATCH -N 1 #number of nodes. Note all your jobs could use a max of 4 nodes on our MBU cluster.
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user ec474@mrc-mbu.cam.ac.uk  #your email

bash /suffolk/Homes/ec474/scripts/count.sh

echo DONE!
