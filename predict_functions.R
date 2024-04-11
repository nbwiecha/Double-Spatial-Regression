# Prediction/estimation functions used in DSR estimation

require(fields)
require(GpGp)
require(mgcv)
source("TV_SVM_function_hpc.R")


spline_fn <- function(Y, S, D, k=100, method="GCV.Cp"){
  trainD.mat <- as.matrix(D)
  p <- ncol(trainD.mat)
  
  fit.Y <- gam(Y ~ D + s(S[,1], S[,2], bs="tp", k=k), method=method)
  
  betahat <- fit.Y$coefficients[2:(p+1)]
  var_betahat <- fit.Y$Vp[2:(p+1), 2:(p+1)]
  
  return(list(betahat, var_betahat))
}
# 



spline_predict <- function(trainS, trainY, trainD, testS, k=100){
  trainD.mat <- as.matrix(trainD)
  
  p <- ncol(trainD.mat)
  
  s1 <- trainS[,1]
  s2 <- trainS[,2]
  fitY <- gam(trainY ~ s(s1, s2, bs="tp", k=k))
  newdata <- data.frame(s1=testS[,1], s2=testS[,2])
  predY <- predict(fitY, newdata=newdata)
  
  predD <- matrix(nrow=nrow(testS), ncol=p)
  for(j in 1:p){
    fit.D.j <- gam(trainD.mat[,j] ~ s(s1, s2, bs="tp", k=k))
    predD[,j] <- predict(fit.D.j, newdata=newdata)
  }
  
  return(list(predY=predY, predD=predD))
}


spatial_plus_fn <- function(Y, D, S, k=100, B=100){
  require(mgcv)
  # Remove spatial trend from treatment to get residuals, then regress Y on residuals
  # Y is response, D treatment, S locations
  dat <- cbind(D, Y, S)
  get_spatialp_estimate <- function(dat, ind, k=k){
    dat <- dat[ind,]
    D <- dat[,1]
    Y <- dat[,2]
    S <- dat[,3:4]
    
    trainD.mat <- as.matrix(D)
    p <- ncol(trainD.mat)
    Rd <- matrix(nrow=nrow(trainD.mat), ncol=p)
    for(j in 1:p){
      fit.D.j <- gam(trainD.mat[,j] ~  s(S[,1], S[,2],bs="tp", k=k))
      Rd[,j] <- residuals(fit.D.j)
    }
    
    fit.Y <- gam(Y ~ Rd + s(S[,1], S[,2], bs="tp", k=k))
    betahat <- fit.Y$coefficients[2]
    return(betahat)
    
  }
  out <- boot(dat, get_spatialp_estimate, R=B, k=k)
  betahat <- out$t0
  var_betahat <- var(out$t)
  
  # se_betahat <- summary(fit.Y)$se[2] # not really correct but whatever
  return(list(betahat, var_betahat))
}

gSEM_fn <- function(Y, D, S, k=100, B=100, include_dsr=FALSE){
  require(mgcv)
  require(boot)
  dat <- cbind(D, Y, S)
  
  get_gsem_estimate <- function(dat, ind, k, include_dsr){
    dat <- dat[ind,]
    
    D <- dat[,1]
    Y <- dat[,2]
    S <- dat[,3:4]
    
    trainD.mat <- as.matrix(D)
    p <- ncol(trainD.mat)
    Rd <- matrix(nrow=nrow(trainD.mat), ncol=p)
    for(j in 1:p){
      fit.D.j <- gam(trainD.mat[,j] ~  s(S[,1], S[,2],bs="tp", k=k))
      Rd[,j] <- residuals(fit.D.j)
    }
    
    fit.Y <- gam(Y ~ s(S[,1], S[,2], bs="tp", k=k))
    Ry <- residuals(fit.Y)
    
    fit.R <- lm(Ry ~ Rd)
    betahat.gsem <- fit.R$coefficients[2]
    
    fit.R.dsr <- lm(Ry ~ Rd - 1)
    betahat.dsr <- fit.R.dsr$coefficients
    if(include_dsr){
      return(c(betahat.gsem, betahat.dsr))
    }else{
      return(betahat.gsem)
    }
    
  }
  
  
  out <- boot(dat, get_gsem_estimate, R=B, k=k, include_dsr=include_dsr)
  
  if(include_dsr){
    betahat.gsem <- out$t0[1]
    var_betahat.gsem <- var(out$t[,1])
    
    betahat.dsr <- out$t0[2]
    var_betahat.dsr <- var(out$t[,2])
    
    return(list(betahat.gsem, var_betahat.gsem, 
                betahat.dsr, var_betahat.dsr))
  }else{
    betahat <- out$t0
    var_betahat <- var(out$t)
    return(list(betahat, var_betahat))
  }
   
}

svm_predict <- function(trainS, trainY, trainD, testS){
  # Prediction function calling training-validation GP estimation function
  n <- nrow(trainS)
  trainD.mat <- as.matrix(trainD)
  p <- ncol(trainD.mat)
  epsilon <- 1/n
  delta <- 1/n^.25
  if(n < 2000){
    delta <- 1/n^.5
  }
  predY <- cv_gp_predict(trainS, testS, trainY, n=n, 
                         epsilon=epsilon, delta=delta, silent=TRUE)$yhat
  predD <- matrix(nrow=nrow(testS), ncol=p)
  for(j in 1:p){
    predD[,j] <- cv_gp_predict(trainS, testS, trainD.mat[,j], n=n,
                                 epsilon=epsilon, delta=delta, silent=TRUE)$yhat
  }
  
  return(list(predY=predY, predD=predD))
}

