#!/bin/tcsh
#BSUB -n 16
#BSUB -W 18:00
#BSUB -R "span[hosts=1]"
#BSUB -J dml_U_gneiting
#BSUB -o stdout.%J
#BSUB -e stderr.%J
cd /share/$GROUP/$USER/DSR
module load R
conda activate env_R421
Rscript run_U_gneiting_sims_hpc.R