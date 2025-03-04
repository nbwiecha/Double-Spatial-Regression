# Prediction/estimation functions used in simulations for Two-Stage Estimators for Spatial Confounding with Point-Referenced Data
# Code by Nate Wiecha, North Carolina State University

require(fields)
require(GpGp)
require(mgcv)
require(randomForest)
library(dbarts)
# setwd("~/GitHub/Spatial-DML")
source("TV_SVM_function_hpc.R")


spatial_reg_gpgp <- function(S,Y,X=NULL,
                             Sp=NULL,
                             Xp=NULL){

  # Spatial regression using GpGp package  
  require(GpGp)
  require(fields)
  n <- length(Y)
  mle <- fit_model(y=Y, locs=S, X=cbind(rep(1, n), X), silent=TRUE)

  yp <- NULL
  if(!is.null(Sp)){
    # If trained on D but no prediction D provided, we want ghat, aka spatial RE's. So create Xp as (1 0)
    if( !is.null(X) & is.null(Xp) ){Xp <- rep(0, nrow(Sp))}
    np <- nrow(Sp)
    yp <- predictions(mle, locs_pred=Sp, X_pred=cbind(rep(1, np),Xp))
  }
  out <- list(mle=mle,pred=yp)
  return(out)}

ols_fn <- function(Y, D){
  # function to get theta hat and se from OLS
  trainD.mat <- as.matrix(D)
  p <- ncol(trainD.mat)
  
  fit <- lm(Y ~ D)
  thetahat <- coefficients(fit)[2:(p+1)]
  var_thetahat <- vcov(fit)[2:(p+1), 2:(p+1)]
  return(list(thetahat, var_thetahat))
}

lmm_fn <- function(Y, S, D){
  # function to get theta hat and se from spatial LMM (using vecchia/GpGp)
  trainD.mat <- as.matrix(D)
  p <- ncol(trainD.mat)
  
  mle <- spatial_reg_gpgp(Y=Y, S=S, X=D)$mle
  thetahat <- mle$betahat[2:(p+1)]
  var_thetahat <- mle$betacov[2:(p+1), 2:(p+1)]
  return(list(thetahat, var_thetahat))
}

spline_fn <- function(Y, S, D, k=100, method="GCV.Cp"){
  trainD.mat <- as.matrix(D)
  p <- ncol(trainD.mat)
  
  fit.Y <- gam(Y ~ D + s(S[,1], S[,2], bs="tp", k=k), method=method)
  
  betahat <- fit.Y$coefficients[2:(p+1)]
  var_betahat <- fit.Y$Vp[2:(p+1), 2:(p+1)]
  
  return(list(betahat, var_betahat))
}


krige_predict_lhat <- function(trainS, trainY, trainD, testS){
  # Kriging prediction function using GpGp for theoretical-style DSR
  # (Which does not include treatment variable in prediction model for outcome)
  trainD.mat <- as.matrix(trainD)
  p <- ncol(trainD.mat)
  
  fit.Y <- spatial_reg_gpgp(S=trainS, Y=trainY, Sp=testS) 
  predY <- fit.Y$pred
  
  predD <- matrix(nrow=nrow(testS), ncol=p)
  for(j in 1:p){
    fit.D.j <- spatial_reg_gpgp(S=trainS, Y=trainD.mat[,j], Sp=testS)
    predD[,j] <- fit.D.j$pred
  }
  
  return(list(predY=predY, predD=predD))
}

spline_predict <- function(trainS, trainY, trainD, testS, k=100){
  
  # Prediction using splines for DSR
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


rf_predict <- function(trainS, trainY, trainD, testS){
  # Prediction using random forest for DSR (not used)
  require(randomForest)
  trainD.mat <- as.matrix(trainD)
  p <- ncol(trainD.mat)
  fit.Y <- randomForest(x=trainS, y=trainY)
  predY <- predict(fit.Y, testS)
  
  predD <- matrix(nrow=nrow(testS), ncol=p)
  for(j in 1:p){
    D.fit <- as.factor(trainD.mat[,j])
    fit.D.j <- randomForest(x=trainS, y=D.fit)
    predD[,j] <- predict(fit.D.j, testS)
  }

  return(list(predY=predY, predD=predD))
}



spatial_plus_fn <- function(Y, D, S, k=100, B=100){
  # Spatial+ estimation
  require(mgcv)
  # Remove spatial trend from treatment to get residuals, then regress Y on residuals
  # Y is response, D treatment, S locations
  dat <- cbind(Y, S, D)
  # p <- ncol(D)
  get_spatialp_estimate <- function(dat, ind, k=k){
    dat <- dat[ind,]
    Y <- dat[,1]
    S <- dat[,2:3]
    D <- dat[,4:ncol(dat)]
    
    trainD.mat <- as.matrix(D)
    p <- ncol(trainD.mat)
    Rd <- matrix(nrow=nrow(trainD.mat), ncol=p)
    for(j in 1:p){
      fit.D.j <- gam(trainD.mat[,j] ~  s(S[,1], S[,2],bs="tp", k=k))
      Rd[,j] <- residuals(fit.D.j)
    }
    
    fit.Y <- gam(Y ~ Rd + s(S[,1], S[,2], bs="tp", k=k))
    betahat <- fit.Y$coefficients[2:(p+1)]
    return(betahat)
    
  }
  out <- boot(dat, get_spatialp_estimate, R=B, k=k)
  betahat <- out$t0
  var_betahat <- var(out$t)
  
  return(list(betahat, var_betahat))
}

gSEM_fn <- function(Y, D, S, k=100, B=100, include_dsr=FALSE){
  # gSEM estimation
  require(mgcv)
  require(boot)
  # dat <- cbind(D, Y, S)
  dat <- cbind(Y, S, D)
  
  get_gsem_estimate <- function(dat, ind, k, include_dsr){
    dat <- dat[ind,]
    
    Y <- dat[,1]
    S <- dat[,2:3]
    D <- dat[,4:ncol(dat)]
    p <- ncol(D)
    
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
    betahat.gsem <- fit.R$coefficients[2:(p+1)]
    
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
    betahat.gsem <- out$t0[1:p]
    var_betahat.gsem <- var(out$t[,1:p])
    
    betahat.dsr <- out$t0[(p+1):(2*p)]
    var_betahat.dsr <- var(out$t[,(p+1):(2*p)])
    
    return(list(betahat.gsem, var_betahat.gsem, 
                betahat.dsr, var_betahat.dsr))
  }else{
    betahat <- out$t0
    var_betahat <- var(out$t)
    return(list(betahat, var_betahat))
  }
  
}

gpgp_predict <- function(trainS, trainY, trainD, testS, fitY=NULL, fitD=NULL, type="ghat",
                         covfun_name="matern_isotropic"){
  
  # type indicates the "type" of DML we're doing (corresponding to the nuisance 
  # parameters and score function used from Chernozhukov 2018), and should be equal 
  # to "ghat" for DSR or "lhat" for the theoretical DSR
  
  trainD.mat <- as.matrix(trainD)
  p <- ncol(trainD.mat)
  
  if(is.null(fitY)){
    if(type=="ghat"){
      fitY <- fit_model(y=trainY, locs=trainS, X=cbind(1,trainD), silent=TRUE,
                        covfun_name=covfun_name)
    }
    
    if(type=="lhat"){
      fitY <- fit_model(y=trainY, locs=trainS, silent=TRUE,
                        covfun_name=covfun_name)
    }
    
  }
  
  if(is.null(fitD)){
    fitD <- list(length=p)
    for(j in 1:p){
      fitD[[j]] <- fit_model(y=trainD.mat[,j], locs=trainS, silent=TRUE,
                             covfun_name=covfun_name)
    }
  }
  
  if(type=="ghat"){
    predY <- predictions(fitY, locs_pred=testS, X_pred=cbind(1, rep(0, nrow(testS))))
  }
  if(type=="lhat"){
    predY <- predictions(fitY, locs_pred=testS, X_pred=rep(1, nrow(testS)))
  }
  
  
  predD <- matrix(nrow=nrow(testS), ncol=ncol(trainD.mat))
  for(j in 1:p){
    predD[,j] <- predictions(fitD[[j]], locs_pred=testS, y_obs=trainD.mat[,j], 
                             locs_obs=trainS,X_obs=rep(1, nrow(trainS)),
                             X_pred=rep(1, nrow(testS)))
  }
  
  return(list(predY=predY, predD=predD))
}

gpgp_predict_X <- function(trainS, trainY, trainD, testS, 
                           trainX, testX){
  
  # Incorporating covariates X
  require(GpGp)
  trainD.mat <- as.matrix(trainD)
  p <- ncol(trainD.mat)
  
  fitY <- spatial_reg_gpgp(S=trainS, Y=trainY, X=trainX,
                           Sp=testS, Xp=testX)
  predY <- fitY$pred
  
  predD <- matrix(nrow=nrow(testS), ncol=ncol(trainD.mat))
  for(j in 1:p){
    fitD <- spatial_reg_gpgp(S=trainS, Y=trainD.mat[,j], X=trainX,
                             Sp=testS, Xp=testX)
    predD[,j] <- fitD$pred
  }
  
  return(list(predY=predY, predD=predD))
}

svm_predict <- function(trainS, trainY, trainD, testS){
  
  # Kriging using the training-validation scheme in Eberts and Steinwart (2013)
  # They analyze this estimator through the lens of support vector machines
  # aka, SVM
  
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


shift_dml = function(Y, D, S, shift=1, B){
  # Shift DML from Gilbert (2021)
  
  # Gilbert's implemented shift estimator function
  # modified to include na.rm=TRUE
  # and to do bootstrapping in the function
  
  suppressMessages(require(dbarts))
  suppressMessages(require(mgcv))
  suppressMessages(require(boot))
  dat <- cbind(D, Y, S)

  get_shift_estimate <- function(dat, ind, shift){
    n <- nrow(dat)
    dat <- dat[ind,]
    X = dat[,1]
    Y = dat[,2]
    lat = dat[,3]
    lon = dat[,4]
    
    df1 = data.frame(lat = lat, lon = lon, X = X)
    df2 = data.frame(lat = lat, lon = lon, X = X+shift)
    mumod = bart(data.frame(lat = lat, lon = lon, X), y.train = Y, x.test = rbind(df2, df1))
    
    mu.fitted.shift = mumod$yhat.test.mean[1:n]
    mu.fitted = mumod$yhat.test.mean[-(1:n)]
    mu_qr = mean(mu.fitted.shift)
    
    mu = mean(Y)
    
    pimod = gam(X ~ s(lat, lon, k=200))
    resid = X - pimod$fitted.values
    d_mod = density(resid)
    
    dens_est = function(x, lat, lon) approx(d_mod$x, d_mod$y, x - predict(pimod, data.frame(lat = lat, lon = lon)))$y
    
    lambda = function(x, lat, lon) dens_est(x-shift, lat, lon)/dens_est(x, lat, lon)
    mu_qipw = sum(lambda(X, lat, lon)*Y, na.rm=TRUE)/sum(lambda(X, lat, lon), na.rm=TRUE)
    
    comp_model = lm(Y ~ lambda(X,lat, lon) + offset(mu.fitted) + 0)
    gamma_hat = coef(comp_model)
    
    mu.hat <- mean(mu.fitted.shift + gamma_hat*lambda(X+shift, lat, lon), na.rm=TRUE) - mu
    
    return(mu.hat)
  }
  
  out <- boot(data=dat, statistic=get_shift_estimate, R=B, shift=shift)
  betahat <- out$t0
  var_betahat <- var(out$t)
  return(list(betahat, var_betahat))
}

