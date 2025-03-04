#!/bin/tcsh
#BSUB -n 32
#BSUB -W 16:00
#BSUB -R "span[hosts=1]"
#BSUB -J dml_deterministic
#BSUB -o stdout.%J
#BSUB -e stderr.%J
cd /share/$GROUP/$USER/DSR
module load R
conda activate env_R421
Rscript run_deterministic_sims_hpc.R