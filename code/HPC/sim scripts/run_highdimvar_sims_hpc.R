# Simulations for high dimensional spatial surface (not used) and higher variance in treatment variable
# For Two-Stage Estimators for Spatial Confounding
# Code by Nate Wiecha
rm(list=ls())

source("dml_function_vector_hpc.R")
source("dml_simulation_function_foreach_hpc.R")
source("dml_simulation_function_foreach_highdim_hpc.R")
source("simulation_functions_hpc.R")


nCores <- strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
cl <- makeCluster(nCores)
registerDoSNOW(cl)

nsims <- 400

# highdim_sim <- run_sim_high_dim(n=1000, nsims=nsims)
highvar_sim <- run_sim_higher_var(n=1000, nsims=nsims)

save(#highdim_sim,
  highvar_sim, file="outputs//sims_highdimvar_n1000.Rdata")
