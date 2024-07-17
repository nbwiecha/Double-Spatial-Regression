# Simulations using a rough (Matern 1.5 covariance) unobserved confounder
# For Two-Stage Estimators for Spatial Confounding
# Code by Nate Wiecha, North Carolina State University
rm(list=ls())

source("dml_function_vector_hpc.R")
source("dml_simulation_function_foreach_hpc.R")
source("simulation_functions_hpc.R")

nCores <- strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
cl <- makeCluster(nCores)
registerDoSNOW(cl)

nsims <- 400

U_matern_sims <- run_sims_U_matern(n=1000, nsims=nsims)

save(U_matern_sims, file="outputs//sims_U_matern_n1000.Rdata")
