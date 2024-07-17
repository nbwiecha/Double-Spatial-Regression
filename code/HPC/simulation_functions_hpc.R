# Full simulation functions for Two-Stage Estimators for Spatial Confounding
# Code by Nate Wiecha, North Carolina State University

get_outputs <- function(result, theta, n){
  bias_thetahat <- mean(result$thetahat - theta)
  rel_bias <- bias_thetahat / theta
  mc_se <- sd(result$thetahat)/sqrt(nrow(result))
  mse_thetahat <- mean((result$thetahat - theta)^2)
  ci_lb <- result$thetahat - 1.96*(result$se_thetahat)
  ci_ub <- result$thetahat + 1.96*(result$se_thetahat)
  ci_lengths <- ci_ub - ci_lb
  mean_length <- mean(ci_lengths)
  contains_theta <- ci_lb < theta & ci_ub > theta
  coverage <- mean(contains_theta)
  power <- mean ( ci_lb > 0 | ci_ub < 0)
  return(c(bias_thetahat, rel_bias, mc_se,mse_thetahat, mean_length, coverage, power))
}

cubed <- function(S, u){return(u^3)}

deterministic <- function(S, u){
  cos(S[,1]*10)*sin(S[,2]*10)
  
}

deterministic2 <- function(S){
  sin(S[,1]*10)*sin(S[,2]*10)
}


gammad <- function(S, u){
  qgamma(pnorm(u/sqrt(3)), 1, 1/sqrt(3))
}

zeroed <- function(S, u){return(0)}

eastwest_error <- function(S, e){
  return(S[,1]*e)
}

middleout_spatial <- function(S, u){
  omega <- pnorm((S[,1]-.5)/.1)
  return(sqrt(omega/3)*u)
}

middleout_error <- function(S, e){
  omega <- pnorm((S[,1]-.5)/.1)
  weight <- sqrt(1-omega)
  return(weight*e)
}

create_output <- function(sim){
  methods <- sim$methods
  output <- matrix(nrow=7, ncol=length(methods))
  rownames(output) <- c("Bias", "Relative Bias", "Monte Carlo SE", "MSE", 
                        "95% CI Length", "95% CI Coverage", "Power")
  colnames(output) <- methods
  
  for(i in 1:length(methods)){
    output[,i] <- get_outputs(sim$results[[i]], sim$theta[1], sim$n)
  }
  
  return(list(
    n=sim$n, nu=sim$nu, confounding_fn=sim$confounding_fn, error_fn=sim$error_fn, output=output, corrs=summary(sim$corrs)
  ))
  
}

run_sim_higherp <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation_p(n=n, nsims, phiD=.2, phiU=.2,
                               nuD=1.5, nuU=1.5, rho=.5, sdY=1, sdD=.1,
                               confounding_fn=NULL, 
                               cov.model.d="matern", cov.model.u="matern",
                               binary_D=binary_D, grid=FALSE,
                               p=3,theta=rep(.5, 3))
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
}

run_sim_deterministic <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation_determ(n=n, nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=0.5, theta=0.5, sdY=1, sdD=.1,
                             confounding_fn=deterministic, y_fn=deterministic2,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D, grid=FALSE)
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
  
}

run_sim_grid <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                                     nuD=NULL, nuU=NULL, rho=0.5, theta=.5, sdY=1, sdD=.1,
                                     confounding_fn=NULL,
                                     cov.model.d="gneiting", cov.model.u="gneiting",
                                     binary_D=binary_D,
                                     grid=TRUE)
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
}

run_sim_nospatial <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                                     nuD=NULL, nuU=NULL, rho=0.5, theta=.5, sdY=1, sdD=.1,
                                     confounding_fn=zeroed,
                                     cov.model.d="gneiting", cov.model.u="gneiting",
                                     binary_D=binary_D)
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
}


run_sim_exponential <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                                     nuD=0.5, nuU=0.5, rho=0.5, theta=.5, sdY=1, sdD=.1,
                                     confounding_fn=NULL,
                                     cov.model.d="matern", cov.model.u="matern",
                                     binary_D=binary_D)
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
}


run_sim_high_dim <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation_highdim(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=0.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=NULL, dim=10, true_dim=2,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D)
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
}

run_sim_higher_var <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=.5, theta=.5, sdY=1, sdD=1,
                             confounding_fn=NULL, dim=2,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D)
  out1 <- create_output(sim1)
  return(list(sim1=sim1, out1=out1))
}

run_sims_unconfounded <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  
  # matern/matern
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=1.5, nuU=1.5, rho=0, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=NULL,
                             cov.model.d="matern", cov.model.u="matern",
                             binary_D=binary_D
  )
  out1 <- create_output(sim1)
  
  # gneiting/gneiting
  sim2 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=0, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=NULL,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D
  )
  out2 <- create_output(sim2)
  
  return(list(sim1=sim1, sim2=sim2,
              out1=out1, out2=out2))
  
}

run_sims_U_matern <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=1.5, nuU=1.5, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=NULL,
                             cov.model.d="matern", cov.model.u="matern",
                             binary_D=binary_D
  )
  out1 <- create_output(sim1)
  
  sim2 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=1.5, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=NULL,
                             cov.model.d="gneiting", cov.model.u="matern",
                             binary_D=binary_D
  )
  out2 <- create_output(sim2)
  
  return(list(sim1=sim1, sim2=sim2,
         out1=out1, out2=out2))
  
}

run_sims_U_gneiting <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=1.5, nuU=NULL, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=NULL,
                             cov.model.d="matern", cov.model.u="gneiting",
                             binary_D=binary_D
  )
  out1 <- create_output(sim1)
  
  sim2 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=NULL,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D
  )
  out2 <- create_output(sim2)
  
  return(list(sim1=sim1, sim2=sim2,
         out1=out1, out2=out2))
}

run_sims_trans_gneiting <- function(n, nsims, error_fn=NULL, binary_D=FALSE){
  sim1 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=1.5, nuU=NULL, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=middleout_spatial,
                             error_fn=middleout_error,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D
  )
  out1 <- create_output(sim1)
  
  sim2 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=cubed,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D
  )
  out2 <- create_output(sim2)
  
  sim3 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=gammad,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D
  )
  
  
  out3 <- create_output(sim3)
  
  sim4 <- run_dml_simulation(n=n, nsims=nsims, phiD=.2, phiU=.2,
                             nuD=NULL, nuU=NULL, rho=.5, theta=.5, sdY=1, sdD=.1,
                             confounding_fn=eastwest_error,
                             cov.model.d="gneiting", cov.model.u="gneiting",
                             binary_D=binary_D
  )
  out4 <- create_output(sim4)
  
  return(list(sim1=sim1, sim2=sim2, sim3=sim3, sim4=sim4,
              out1=out1, out2=out2, out3=out3, out4=out4))
}
