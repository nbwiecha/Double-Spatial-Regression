#!/bin/tcsh
#BSUB -n 16
#BSUB -W 18:00
#BSUB -R "span[hosts=1]"
#BSUB -J dml_highdimvar
#BSUB -o stdout.%J
#BSUB -e stderr.%J
cd /share/$GROUP/$USER/DSR
module load R
conda activate env_R421
Rscript run_highdimvar_sims_hpc.R