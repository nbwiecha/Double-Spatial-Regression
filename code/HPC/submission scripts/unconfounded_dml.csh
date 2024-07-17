#!/bin/tcsh
#BSUB -n 16
#BSUB -W 18:00
#BSUB -R "span[hosts=1]"
#BSUB -J dml_unconfounded
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -q stat
cd /share/$GROUP/$USER/DML_sims
module load R
conda activate env_R421
Rscript run_unconfounded_sims_hpc.R