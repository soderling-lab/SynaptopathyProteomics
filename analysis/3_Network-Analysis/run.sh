#!/bin/bash
#SBATCH -a 1          # How many times you execute? e.g. 1-4
#SBATCH -e out.stderr # Where is stderr sent?
#SBATCH -o out.stdout # Where is stdout sent?
#SBATCH -N 1          # How many nodes?
#SBATCH -c 16          # How many cores?
#SBATCH --mem=16G      # How much memory to allocate?
#SBATCH -p scavenger  # scavenge some CPU.
module load R/3.6.0
Rscript 4_self-preservation.R
#R CMD BATCH test.R


