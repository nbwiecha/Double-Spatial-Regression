# Two-Stage Estimators for Spatial Confounding with Point-Referenced Data
# Code by Nate Wiecha, North Carolina State University

# Run simulations: scenarios where latent functions of space are determined deterministically.

rm(list=ls())

# setwd("~/GitHub/Spatial-DML/code/HPC")
source("dml_function_vector_hpc.R")
source("dml_simulation_function_foreach_determ_hpc.R")
# source("dml_simulation_function_foreach_highdim_hpc.R")
source("simulation_functions_hpc.R")

nCores <- strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
cl <- makeCluster(nCores)
# cl <- makeCluster(7)
registerDoSNOW(cl)

nsims <- 400

run_sim_deterministic1 <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  # Alternate function where y_fn=NULL so same fn of space in D and Y
  sim1 <- run_dml_simulation_determ(n=n, nsims, phiD=.2, phiU=.2,
                                    nuD=NULL, nuU=NULL, rho=0.5, theta=0.5, sdY=1, sdD=.1,
                                    confounding_fn=deterministic, y_fn=NULL,
                                    cov.model.d="gneiting", cov.model.u="gneiting",
                                    binary_D=binary_D, grid=FALSE)
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
  
}
deterministic_sim_diff <- run_sim_deterministic(n=1000, nsims=nsims)
deterministic_sim_same <- run_sim_deterministic1(n=1000, nsims=nsims)
save(deterministic_sim_same, deterministic_sim_diff, file="outputs//sims_deterministic_n1000.Rdata")

