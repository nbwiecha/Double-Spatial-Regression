# Two-Stage Estimators for Spatial Confounding with Point-Referenced Data
# Code by Nate Wiecha, North Carolina State University

# Run simulations: scenarios where latent functions of space are transformed.
.
rm(list=ls())

# setwd("~/GitHub/Spatial-DML/code/HPC")
source("dml_function_vector_hpc.R")
source("dml_simulation_function_foreach_hpc.R")
source("simulation_functions_hpc.R")

nCores <- strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
cl <- makeCluster(nCores)
registerDoSNOW(cl)

nsims <- 400

trans_sims <- run_sims_trans_gneiting(n=1000, nsims=nsims)

save(trans_sims, file="outputs//sims_trans_n1000.Rdata")
