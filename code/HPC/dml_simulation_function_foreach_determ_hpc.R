# Function to perform simulation study for Two-Stage Estimators for Spatial Confounding with Point-Referenced Data
# Simulations for deterministic spatial surfaces in D, Y as opposed to random GP draws
# Code by Nate Wiecha, North Carolina State University

library(foreach)
library(doSNOW)
library(doParallel)
library(fields)
library(geoR)
library(boot)
# setwd("~/GitHub/Spatial-DML")

run_dml_simulation_determ <- function(
                               n, nsims,     # sample size, length of D, number of MC samples
                               grid=FALSE,   # Whether to use regular grid of points
                               randomfields=FALSE, # whether to use random fields to generate data
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
                               y_fn=NULL, # function of space to add arbitrary function of space into Y
                               error_fn=NULL,#, # function of iid error that goes into Y, eg nonstationary
                               u.df=Inf, # df for u's t-distn (Inf=>normal)
                               K=5, # number of folds to use for DML cross-fitting
                               B=100, # number of bootstrap replicates for gsem and spatial+
                               dim=2 # dimension of spatial domain; only for grid=FALSE
){
  require(randomForest)
  require(fields)
  require(GpGp)
  require(RandomFields)
  require(mgcv)
  require(geoR)
  
  p <- 1 # dimension of treatment
  
  # Set up MC simulation
  methods <- c("OLS", "Spline.GCV", "Spline.REML", "LMM", "SpatialPlus", 
               # "gSEM_int", 
               "gSEM_noint",
               "Shift_DML",  
               "DML_GP", "DML_GpGp",  
               "DML_alt", 
               "DML_GP_nosplit", "DML_alt_nosplit",
               "DML_spline_nosplit", "DML_spline_nosplit_lhat",
               "DML_lhat_smooth", "DML_alt_smooth")

  
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
  run_sim_iter <- function(methods, seed){
    set.seed(seed)
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
    
    if(is.null(y_fn)){
      y_fn2 <- function(S){0}
    }else{
      y_fn2 <- y_fn
    }
    source("dml_function_vector_hpc.R", local=TRUE)
    source("dml_function_vector_nosplit_hpc.R", local=TRUE)
    source("predict_functions_hpc.R", local=TRUE)
    
    iter_result <- matrix(nrow=length(methods), ncol=2)
    rownames(iter_result) <- methods
    colnames(iter_result) <- c("thetahat", 'se_thetahat')
   
    

    # Old code using RandomFields: will be needed for large n
    
    if(randomfields){
      x <- seq(0, 1, length=sqrt(n))
      y <- seq(0, 1, length=sqrt(n))
      S <- as.matrix(expand.grid(x,y))

      # Set up Matern model for D, U
      rho_mat <- diag(p+1)
      rho_mat[(1:p),(p+1)] <- rho
      rho_mat[(p+1), 1:p] <- rho
      
      nu <- c(nuD, nuU)
      #
      model.DU <- RMparswmX(nudiag=nu, rho=rho_mat, scale=phiDU)
      #
      DU <- RFsimulate(model.DU, x=x, y=y)
      D <- as.matrix(DU@data[,1:p])
      for(j in 1:p){
        D[,j] <- D[,j] + rnorm(n=n, sd=sdD)
      }
      U <- as.matrix(DU@data[,(p+1)])
      
      if (u.df!=Inf){
        U <- U / sqrt(rchisq(n, u.df)/u.df)
      }
    }else{
      if(grid){
        vecs <- list(length=dim)
        for(dim.i in 1:dim){
          vecs[[dim.i]] <- seq(0, 1, length=n^{1/dim})
        }
        
        S <- as.matrix(
          expand.grid(
            vecs
          )
        )
      }else{
        S <- matrix(data=runif(n * dim), nrow=n)
      }
      
      d <- rdist(S) 
      R.D <- cov.spatial(d, cov.model=cov.model.d, cov.pars=c(1, phiD), kappa=nuD)
      R.U <- cov.spatial(d, cov.model=cov.model.u, cov.pars=c(1, phiU), kappa=nuU)
      # dim(R.D)
      
      L.D <- chol(R.D)

      L.U <- chol(R.U)

      Sigma1 <- cbind(R.D, rho * t(L.D) %*% L.U)
      Sigma2 <- cbind(rho * t(L.U) %*% L.D, R.U)
      Sigma <- rbind(Sigma1, Sigma2)

      D.U <- rmvnorm(n=1, sigma=Sigma)
      D <- as.matrix(h(S, D.U[1:n]) + rnorm(n, sd=sdD))
      U <- D.U[(n+1):(2*n)]
      if (u.df!=Inf){
        U <- U / sqrt(rchisq(n, u.df)/u.df)
      }
    }
    
    Y <- D * theta + h(S, U) + y_fn2(S) + g(S, rnorm(n=n, sd=sdY))
    D <- D - mean(D)
    Y <- Y - mean(Y)
    

    if(p==1){
      iter_result["OLS",] <- as.numeric(ols_fn(Y=Y, D=D))
      iter_result["Spline.GCV",] <- as.numeric(spline_fn(Y=Y, D=D, S=S, k=300))
      iter_result["Spline.REML",] <- as.numeric(spline_fn(Y=Y, D=D, S=S, method="REML", k=300))
      iter_result["LMM",] <- as.numeric(lmm_fn(Y, S, D))
      iter_result["SpatialPlus",] <- as.numeric(spatial_plus_fn(Y, D, S, k=300, B=B))
      # iter_result["gSEM_k100",] <- as.numeric(gSEM_fn(Y, D, S))
      gsem.fit <- gSEM_fn(Y, D, S, k=300, B=B, include_dsr=TRUE)
      # iter_result["gSEM_int",] <- as.numeric(gsem.fit[1:2])
      iter_result["gSEM_noint",] <- as.numeric(gsem.fit[3:4])
      iter_result["Shift_DML",] <- tryCatch(
        {
          as.numeric(shift_dml(Y, D, S, B=0))
        },
        error=function(e){return(NA)})
      iter_result["DML_GP",] <- as.numeric(dml_lhat(Y, D, S, K=K, predict_fn=svm_predict))
      iter_result["DML_GpGp",] <- as.numeric(dml_lhat(Y, D, S, K=K, predict_fn=gpgp_predict, type="lhat"))
      iter_result["DML_alt",] <- as.numeric(dml_alt(Y, D, S, X=NULL, K=K, refit=TRUE))
      iter_result["DML_GP_nosplit",] <- as.numeric(dml_lhat_nosplit(Y, D, S, predict_fn=svm_predict))
      iter_result["DML_alt_nosplit",] <- as.numeric(dml_alt_nosplit(Y, D, S, X=NULL))
      iter_result["DML_spline_nosplit",] <- as.numeric(dml_alt_nosplit_spline(Y, D, S, k=300))
      iter_result["DML_spline_nosplit_lhat",] <- as.numeric(dml_lhat_nosplit(Y, D, S, 
                                                                            predict_fn=spline_predict, k=300))
      iter_result["DML_lhat_smooth",] <- as.numeric(dml_lhat(Y, D, S, K=K, predict_fn=gpgp_predict, type="lhat",
                                                             covfun_name="matern45_isotropic"))
      iter_result["DML_alt_smooth",] <- as.numeric(dml_alt(Y, D, S, X=NULL, K=K, refit=TRUE,
                                                           covfun_name="matern45_isotropic"))
    }else{
      print("p>1 not implemented yet")
      return()
    }
    iter_result[,2] <- sqrt(iter_result[,2])
    return(iter_result)
  }
  
  results_list1 <- foreach(i=1:nsims,
                           .packages=c('geoR', 'RandomFields','GpGp','mgcv','fields', 'boot', 'dbarts')) %dopar% 
    run_sim_iter(methods, i)
  
  for(i in 1:nsims){
    for(method in methods)
      results[[method]][i,] <- results_list1[[i]][method,]
  }
  
  return(
    list(results=results, n=n, phiDU=phiDU, nuD=nuD, nuU=nuU, corrs=corrs,
         theta=theta, methods=methods, sds=c(sdY, sdD), rho=rho,
         confounding_fn=confounding_fn,
         error_fn=error_fn, cov.model.d=cov.model.d,grid=grid,
         phiD=phiD,phiU=phiU)
  )
  
}

