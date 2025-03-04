# Function to perform DML simulation for easier use
# Nate Wiecha 1/13/2023
# Updated 7/4/2023: Updating to use DML function for vector D; 
# incorporate setting of whether to do binary D or not, truncating inside DML function;
# moved all the helper functions to separate script;

# Update 10/9/2023:
# add the JHU implemented DML shift estimator

# 2/28/2024: 
# add option for higher-dimensional spatial domain

# Currently does not incorporate spatial REs into Y
library(foreach)
library(doSNOW)
library(doParallel)
library(fields)
library(geoR)
library(boot)
# setwd("~/GitHub/Spatial-DML")

run_dml_simulation_highdim <- function(
                               n, nsims,     # sample size, length of D, number of MC samples
                               # grid=FALSE,   # Whether to use regular grid of points
                               phiDU=NULL,        # Spatial range of D, U (same)
                               phiD=0.2,
                               phiU=0.2,
                               nuD, nuU,     # Smoothness parameters for D, U
                               rho,          # corr. param.s bw D, U (not actual corr)
                               theta,        # treatment effect of D on Y; p-vector
                               sdY=1, sdD=1, # nugget variance in Y and D
                               # progress_bar=TRUE, # Progress bar
                               # message=FALSE,
                               cov.model.d="matern", # covariance function if not using matern
                               cov.model.u="gneiting",
                               binary_D=FALSE,
                               confounding_fn=NULL, # function of unobserved confounder that goes into Y
                               error_fn=NULL,#, # function of iid error that goes into Y, eg nonstationary
                               u.df=Inf, # df for u's t-distn (Inf=>normal)
                               K=5, # number of folds to use for DML cross-fitting
                               B=100, # number of bootstrap replicates for gsem and spatial+
                               dim=2, # dimension of spatial domain; only for grid=FALSE
                               true_dim=2 # number of dimensions of spatial domain used for U
                               # fixed_smoothness=NULL # If want to run simulation without estimating smoothness
){
  require(randomForest)
  require(fields)
  require(GpGp)
  require(RandomFields)
  require(mgcv)
  require(geoR)
  
  p <- 1 # dimension of treatment
  
  
  # Set up grid of spatial locations
  # S      <- as.matrix(expand.grid(seq(0, 1, length=sqrt(n)), seq(0, 1, length=sqrt(n))))
  # x <- y <- seq(0, 1, length=sqrt(n))
  
  
  # Set up MC simulation
  methods <- c("OLS",  "LMM",  "DML_alt", "DML_alt_nosplit",
              "DML_SV_nosplit", "DML_SV")
  # if(!LOO){methods <- c("OLS", "Spline", "SpatialPlus", "gSEM",
                        # "DML_GP_ghat", 
                        # "DML_spline_lhat",  "DML_RF")}
  
  ## Store results as list of matrices which contain the MC estimates of 
  ## theta hat and se(theta hat)
  result1 <- matrix(nrow=nsims, ncol=2)
  result <- as.data.frame(result1)
  colnames(result) <- c("thetahat", "se_thetahat")
  
  results <- list(length=length(methods))
  for(i in 1:length(methods)){
    results[[i]] <- result
  }
  names(results) <- methods
  
  # store sample correlations between D, U: since rho is not the actual corr.
  corrs <- rep(NA, nsims) 
  
  # Run the simulation
  # if(progress_bar){pb <- txtProgressBar(min=0, max=nsims, style=3)}
  nu <- c(nuD, nuU)
  run_sim_iter <- function(methods){
    
    if(is.null(confounding_fn)){
      h <- function(s, u){u}
    }else{
      h <- confounding_fn
    }
    
    if(is.null(error_fn)){
      g <- function(S, e){e}
    }else{
      g <- error_fn
    }
    source("dml_function_vector_hpc.R", local=TRUE)
    source("dml_function_vector_nosplit_hpc.R", local=TRUE)
    source("predict_functions_hpc.R", local=TRUE)
    source("vecchia_scaled.R", local=TRUE)
    
    iter_result <- matrix(nrow=length(methods), ncol=2)
    rownames(iter_result) <- methods
    colnames(iter_result) <- c("thetahat", 'se_thetahat')
   
    S <- matrix(data=runif(n * dim), nrow=n)
    S_use <- S[,1:true_dim]
    # x <- runif(n=n)
    # y <- runif(n=n)
    # S <- cbind(x,y)
    
    d <- rdist(S_use) 
    R.D <- cov.spatial(d, cov.model=cov.model.d, cov.pars=c(1, phiD), kappa=nuD)
    R.U <- cov.spatial(d, cov.model=cov.model.u, cov.pars=c(1, phiU), kappa=nuU)
    # dim(R.D)
    
    L.D <- chol(R.D)
    # ev.D <- eigen(R.D, symmetric=TRUE)
    # L.D <- t(ev.D$vectors %*% (t(ev.D$vectors) * sqrt(pmax(ev.D$values, 0))))
    
    L.U <- chol(R.U)
    # ev.U <- eigen(R.U, symmetric=TRUE)
    # L.U <- t(ev.U$vectors %*% (t(ev.U$vectors) * sqrt(pmax(ev.U$values, 0))))
    
    Sigma1 <- cbind(R.D, rho * t(L.D) %*% L.U)
    Sigma2 <- cbind(rho * t(L.U) %*% L.D, R.U)
    Sigma <- rbind(Sigma1, Sigma2)
    # L.Sigma <- chol(Sigma)
    # D.U <- t(L.Sigma) %*% rnorm(2*n)
    D.U <- rmvnorm(n=1, sigma=Sigma)
    D <- as.matrix(D.U[1:n] + rnorm(n, sd=sdD))
    U <- D.U[(n+1):(2*n)]
    if (u.df!=Inf){
      U <- U / sqrt(rchisq(n, u.df)/u.df)
    }
    
    Y <- D * theta + h(S, U) + g(S, rnorm(n=n, sd=sdY))
    D <- D - mean(D)
    Y <- Y - mean(Y)
    # start_covparms_Y <- get_start_parms(y=Y, X=cbind(rep(1,n)), locs=S, covfun_name="matern_isotropic")
    # start_covparms_D <- numeric(length=p)
    # for(j in 1:p){
    #   start_covparms_D[j] <- get_start_parms(y=D, X=cbind(rep(1,n)), locs=S, covfun_name="matern_isotropic")
    # }
# 
#     fitY <- fit_model(y=Y, locs=S, X=cbind(rep(1, n)), silent=TRUE, #fixed_parms=c(3),
#                       start_parms=start_covparms_Y$start_parms)
# 
#     fitD <- fit_model(y=D, locs=S, X=cbind(rep(1, n)), silent=TRUE, #fixed_parms=c(3),
#                       start_parms=start_covparms_D[j]$start_parms)
#     
    # Below is for vector D; not used in sims
    # fitD <- list(length=p)
    # for(j in 1:p){
    #   fitD[[j]] <- fit_model(y=D[,j], locs=S, X=cbind(rep(1, n)), silent=TRUE, #fixed_parms=c(3),
    #                         start_parms=start_covparms_D[j]$start_parms)
    # }

    if(p==1){
      iter_result["OLS",] <- as.numeric(ols_fn(Y=Y, D=D))
      # iter_result["Spline.GCV",] <- as.numeric(spline_fn(Y=Y, D=D, S=S))
      # iter_result["Spline.REML",] <- as.numeric(spline_fn(Y=Y, D=D, S=S, method="REML"))
      iter_result["LMM",] <- as.numeric(lmm_fn(Y, S, D))
      # iter_result["SpatialPlus",] <- as.numeric(spatial_plus_fn(Y, D, S, k=300, B=B))
      # iter_result["gSEM_k100",] <- as.numeric(gSEM_fn(Y, D, S))
      # gsem.fit <- gSEM_fn(Y, D, S, k=300, B=B, include_dsr=TRUE)
      # iter_result["gSEM_int",] <- as.numeric(gsem.fit[1:2])
      # iter_result["gSEM_noint",] <- as.numeric(gsem.fit[3:4])
      # iter_result["Shift_DML",] <- tryCatch(
      #   {
      #     as.numeric(shift_dml(Y, D, S, B=0))
      #   },
      #   error=function(e){return(NA)})
      # iter_result["DML_spline_k100",] <- as.numeric(dml_lhat(Y, D, S, K=5, predict_fn=spline_predict))
      # iter_result["DML_spline",] <- as.numeric(dml_lhat(Y, D, S, K=5, predict_fn=spline_predict, k=300))
      # iter_result["DML_GP",] <- as.numeric(dml_lhat(Y, D, S, K=K, predict_fn=svm_predict))
      # iter_result["DML_RF",] <- as.numeric(dml_lhat(Y, D, S, K=5, predict_fn=rf_predict))
      # iter_result["DML_GpGp",] <- as.numeric(dml_lhat(Y, D, S, K=5, predict_fn=gpgp_predict,
      #                                                 fitY=fitY, fitD=fitD))
      
      # The scaled vecchia function throws an overflow error "randomly"
      highdim.fit <- tryCatch(
        {
          dml_alt_highdim(Y, D, S, K=K)
        },
        error=function(e){
          message("Trying again")
          tryCatch(
            {
              dml_alt_highdim(Y, D, S, K=K)
            },
            error=function(e){
              message("no")
              return(c(NA,NA))
            }
          )
          
        }
      )
      
      highdim.fit.nosplit <- tryCatch(
        {
          dml_alt_highdim_nosplit(Y, D, S)
        },
        error=function(e){
          message("Trying again")
          tryCatch({
            dml_alt_highdim_nosplit(Y, D, S)
          },
          error=function(e){
            message("no")
            return(c(NA,NA))
          }
          )
        }
      )
      iter_result["DML_alt",] <- as.numeric(dml_alt(Y, D, S, X=NULL, K=K))
      iter_result["DML_alt_nosplit",] <- as.numeric(dml_alt_nosplit(Y, D, S, X=NULL))
      
      iter_result["DML_SV",] <- as.numeric(highdim.fit)
      # iter_result["DML_GP_nosplit",] <- as.numeric(dml_lhat_nosplit(Y, D, S, predict_fn=svm_predict))
      iter_result["DML_SV_nosplit",] <- as.numeric(highdim.fit.nosplit)
      # iter_result["DML_spline_nosplit",] <- as.numeric(dml_alt_nosplit_spline(Y, D, S, k=300))
      # iter_result["DML_spline_nosplit_lhat",] <- as.numeric(dml_lhat_nosplit(Y, D, S, 
                                                                            # predict_fn=spline_predict, k=300))
    }else{
      print("p>1 not implemented yet")
      return()
    }
    iter_result[,2] <- sqrt(iter_result[,2])
    return(iter_result)
  }
  
  results_list1 <- foreach(i=1:nsims,
                           .packages=c('geoR', 'RandomFields','GpGp','mgcv','fields', 'boot', 'dbarts')) %dopar% 
    run_sim_iter(methods)
  
  for(i in 1:nsims){
    for(method in methods)
      results[[method]][i,] <- results_list1[[i]][method,]
      # results[[method]][i,2] <- sqrt(results[[method]][i,2]) # double check at some point
  }
  
  return(
    list(results=results, n=n, phiDU=phiDU, nuD=nuD, nuU=nuU, corrs=corrs,
         theta=theta, methods=methods, sds=c(sdY, sdD), rho=rho,
         confounding_fn=confounding_fn,
         error_fn=error_fn, cov.model.d=cov.model.d,grid=grid,
         phiD=phiD,phiU=phiU)
  )
  
}

