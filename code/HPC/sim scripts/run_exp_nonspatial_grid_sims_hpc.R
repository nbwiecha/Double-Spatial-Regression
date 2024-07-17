# Run simulations with exponential covariance, no spatial variation, gridded locations
# For Two-Stage Estimators for Spatial Confounding
# Code by Nate Wiecha, North Carolina State University
rm(list=ls())

source("dml_function_vector_hpc.R")
source("dml_simulation_function_foreach_hpc.R")
source("dml_simulation_function_foreach_highdim_hpc.R")
source("simulation_functions_hpc.R")


nCores <- strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
cl <- makeCluster(nCores)
registerDoSNOW(cl)

set.seed(1234)
nsims <- 400

exponential_sim <- run_sim_exponential(n=1000, nsims=nsims)
nonspatial_sim <- run_sim_nospatial(n=1000, nsims=nsims)
grid_sim <- run_sim_grid(n=1024, nsims=nsims)

save(exponential_sim,nonspatial_sim, grid_sim, file="outputs//sims_exp_nonspatial_grid_n1000.Rdata")
