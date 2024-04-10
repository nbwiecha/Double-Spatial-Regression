# Generic DML Function for Partially Linear Model without cross-fitting
# Nate Wiecha

# predict_fn will need to be specified, it makes predictions for Y and D
# Y is response, D treatment (vector), S spatial locations, X covariate
# K is number of folds, folds is vector of fold assignments

# predict_fn will take as input trainY, trainD, trainS, testS, (trainX, testX)
# and output a list with components predY, predD

# _lhat corresponds to score (4.4) from Chernozhukov
# _alt corresponds to score (4.3) from Chernozhukov

# See also: DoubleML R package which is the official implementation of Chernozhukov (2018)

dml_lhat_nosplit <- function(Y, D, S, predict_fn, ...){
  
  truncate <-  function(x, a, b) {
    x[x < a] <- a
    x[x > b] <- b
    x
  }
  
  psi_a <- function(D, Dhat){
    n <- nrow(D)
    p <- ncol(D)
    arr <- array(dim=c(n, p, p))
    for(i in 1:n){
      arr[i,,] <- -(D - Dhat)[i,] %*% t((D - Dhat)[i,])
    }
    return(arr)
  }
  
  psi_b <- function(Y, D, Yhat, Dhat){
    (D - Dhat) * as.numeric(Y - Yhat)
  }
  
  psi <- function(Y, D, Yhat, Dhat, thetahat){
    
    psi.a <- psi_a(D, Dhat)
    psi.b <- psi_b(Y, D, Yhat, Dhat)

    n <- nrow(psi.b)
    psi <- matrix(nrow=nrow(psi.b), ncol=ncol(psi.b))
        
    for(i in 1:n){
      psi[i,] <- psi.a[i,,] %*% thetahat + psi.b[i,]
    }
    return(psi)
  }
  
  n <- length(Y)
  p <- ncol(D)

  # Put the predictions for each fold in lists
  predictions <- predict_fn(S, Y, D, S, ...)
  predictionsY <- as.numeric(predictions$predY)
  predictionsD <- as.numeric(predictions$predD)
  for(j in 1:p){
    if( (min(as.matrix(D)[,j]) == 0) & (max(as.matrix(D)[,j])== 1) ){
      predictionsD[,j] <- truncate(predictionsD[,j], 0, 1)
    }
  }
  
  # Put the list of predictions into vectors
  Yhat <- predictionsY
  Dhat <- predictionsD
  
  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- lm((Y-Yhat) ~ (Vhat)-1)$coefficients#solve(t(Vhat) %*% (Vhat)) %*% (t(Vhat) %*% (Y - Yhat))
  # Get variance estimate
  Jhat <- colMeans(psi_a(D, Dhat))
  psi <- psi(Y, D, Yhat, Dhat, thetahat)
  Bhat <- 1/n*t(psi) %*% psi
  var_thetahat <- solve(Jhat) %*% Bhat %*% t(solve(Jhat)) / n

  return(list(thetahat, var_thetahat))
}


################################################################################
#                Alternative DSR estimator - GpGp only                         #
################################################################################

dml_alt_nosplit <- function(Y, D, S, X){
  require(GpGp)
  truncate <-  function(x, a, b) {
    x[x < a] <- a
    x[x > b] <- b
    x
  }
  
  psi_a <- function(D, Dhat){
    n <- nrow(D)
    p <- ncol(D)
    arr <- array(dim=c(n, p, p))
    for(i in 1:n){
      arr[i,,] <- -(D - Dhat)[i,] %*% t((D)[i,])
    }
    return(arr)
  }
  
  psi_b <- function(Y, D, Yhat, Dhat){
    (D - Dhat) * as.numeric(Y - Yhat)
  }
  
  psi <- function(Y, D, Yhat, Dhat, thetahat){
    
    psi.a <- psi_a(D, Dhat)
    psi.b <- psi_b(Y, D, Yhat, Dhat)
    
    n <- nrow(psi.b)
    psi <- matrix(nrow=nrow(psi.b), ncol=ncol(psi.b))
    
    for(i in 1:n){
      psi[i,] <- psi.a[i,,] %*% thetahat + psi.b[i,]
    }
    return(psi)
  }
  
  n <- length(Y)
  p <- ncol(D)

  # estimate coefficients on full sample for speed
  fitY <- fit_model(y=Y, locs=S, X=cbind(1, D, X),silent=TRUE)
  fitD <- list(length=p)
  for(j in 1:p){
    fitD[[j]] <- fit_model(y=D[,j], locs=S, X=cbind(rep(1,n),X), silent=TRUE)
  }
  # Put the predictions for each fold in lists
  predictionsY <- predictions(fitY, locs_pred=S,
                                   X_pred=cbind(1, matrix(0, nrow=nrow(S), ncol=p), X))
  predD <- matrix(nrow=nrow(S), ncol=p)
  
  for(j in 1:p){
    predD[,j] <- predictions(fitD[[j]], locs_pred=S,
                               X_pred=cbind(rep(1,n), X))
  }
  predictionsD <- predD
  
  for(j in 1:p){
    if( (min(as.matrix(D)[,j]) == 0) & (max(as.matrix(D)[,j])== 1) ){
      predictionsD[,j] <- truncate(predictionsD[,j], 0, 1)
    }
  }
  
  
  # Put the list of predictions into vectors
  Yhat <- predictionsY
  Dhat <- predictionsD

  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- solve(t(Vhat) %*% (D)) %*% (t(Vhat) %*% (Y - Yhat))
  # Get variance estimate
  Jhat <- colMeans(psi_a(D, Dhat))
  psi <- psi(Y, D, Yhat, Dhat, thetahat)
  Bhat <- 1/n*t(psi) %*% psi
  var_thetahat <- solve(Jhat) %*% Bhat %*% t(solve(Jhat)) / n

  return(list(thetahat, var_thetahat))
}

################################################################################
#                Alternative DSR estimator - spline only                       #
################################################################################

# no covariates for this
dml_alt_nosplit_spline <- function(Y, D, S, k=100){
  require(mgcv)
  truncate <-  function(x, a, b) {
    x[x < a] <- a
    x[x > b] <- b
    x
  }
  
  psi_a <- function(D, Dhat){
    n <- nrow(D)
    p <- ncol(D)
    arr <- array(dim=c(n, p, p))
    for(i in 1:n){
      arr[i,,] <- -(D - Dhat)[i,] %*% t((D)[i,])
    }
    return(arr)
  }
  
  psi_b <- function(Y, D, Yhat, Dhat){
    (D - Dhat) * as.numeric(Y - Yhat)
  }
  
  psi <- function(Y, D, Yhat, Dhat, thetahat){
    
    psi.a <- psi_a(D, Dhat)
    psi.b <- psi_b(Y, D, Yhat, Dhat)
    
    n <- nrow(psi.b)
    psi <- matrix(nrow=nrow(psi.b), ncol=ncol(psi.b))
    
    for(i in 1:n){
      psi[i,] <- psi.a[i,,] %*% thetahat + psi.b[i,]
    }
    return(psi)
  }
  
  n <- length(Y)
  p <- ncol(D)
  
  # estimate coefficients on full sample for speed
  # currently with no X included
  gam.dat <- data.frame(Y, D, S)
  colnames(gam.dat) <- c("Y", "D", "s1", "s2")
  fitY <- gam(Y ~ D + s(s1, s2, k=k), data=gam.dat)
  fitD <- list(length=p)
  for(j in 1:p){
    fitD[[j]] <- gam(D ~ s(s1, s2, k=k), data=gam.dat)
  }
  # Put the predictions for each fold in lists
  pred.gam.dat <- data.frame(Y, D=0, S)
  colnames(pred.gam.dat) <- c("Y", "D", "s1", "s2")
  predictionsY <- matrix(predict(fitY, newdata=pred.gam.dat), nrow=n)
  predD <- matrix(nrow=nrow(S), ncol=p)
  
  for(j in 1:p){
    predD[,j] <- predict(fitD[[j]], newdata=pred.gam.dat)
  }
  predictionsD <- predD
  
  for(j in 1:p){
    if( (min(as.matrix(D)[,j]) == 0) & (max(as.matrix(D)[,j])== 1) ){
      predictionsD[,j] <- truncate(predictionsD[,j], 0, 1)
    }
  }
  
  # Put the list of predictions into vectors
  Yhat <- predictionsY
  Dhat <- predictionsD
  
  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- solve(t(Vhat) %*% (D)) %*% (t(Vhat) %*% (Y - Yhat))
  # Get variance estimate
  Jhat <- colMeans(psi_a(D, Dhat))
  psi <- psi(Y, D, Yhat, Dhat, thetahat)
  Bhat <- 1/n*t(psi) %*% psi
  var_thetahat <- solve(Jhat) %*% Bhat %*% t(solve(Jhat)) / n

  return(list(thetahat, var_thetahat))
}

################################################################################
#                Alternative DSR estimator - GpGp only                         #
################################################################################

dml_alt_nosplit <- function(Y, D, S, X){
  require(GpGp)
  truncate <-  function(x, a, b) {
    x[x < a] <- a
    x[x > b] <- b
    x
  }
  
  psi_a <- function(D, Dhat){
    n <- nrow(D)
    p <- ncol(D)
    arr <- array(dim=c(n, p, p))
    for(i in 1:n){
      arr[i,,] <- -(D - Dhat)[i,] %*% t((D)[i,])
    }
    return(arr)
  }
  
  psi_b <- function(Y, D, Yhat, Dhat){
    (D - Dhat) * as.numeric(Y - Yhat)
  }
  
  psi <- function(Y, D, Yhat, Dhat, thetahat){
    
    psi.a <- psi_a(D, Dhat)
    psi.b <- psi_b(Y, D, Yhat, Dhat)
    
    n <- nrow(psi.b)
    psi <- matrix(nrow=nrow(psi.b), ncol=ncol(psi.b))
    
    for(i in 1:n){
      psi[i,] <- psi.a[i,,] %*% thetahat + psi.b[i,]
    }
    return(psi)
  }
  
  n <- length(Y)
  p <- ncol(D)

  # estimate coefficients and GP params on full sample for speed
  fitY <- fit_model(y=Y, locs=S, X=cbind(1, D, X),silent=TRUE)
  fitD <- list(length=p)
  for(j in 1:p){
    fitD[[j]] <- fit_model(y=D[,j], locs=S, X=cbind(rep(1,n),X), silent=TRUE)
  }
  # Put the predictions for each fold in lists
  predictionsY <- predictions(fitY, locs_pred=S,
                                   X_pred=cbind(1, matrix(0, nrow=nrow(S), ncol=p), X))
  predD <- matrix(nrow=nrow(S), ncol=p)
  
  for(j in 1:p){
    predD[,j] <- predictions(fitD[[j]], locs_pred=S,
                               X_pred=cbind(rep(1,n), X))
  }
  predictionsD <- predD
  
  for(j in 1:p){
    if( (min(as.matrix(D)[,j]) == 0) & (max(as.matrix(D)[,j])== 1) ){
      predictionsD[,j] <- truncate(predictionsD[,j], 0, 1)
    }
  }
  
  
  # Put the list of predictions into vectors
  Yhat <- predictionsY
  Dhat <- predictionsD

  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- solve(t(Vhat) %*% (D)) %*% (t(Vhat) %*% (Y - Yhat))
  # Get variance estimate
  Jhat <- colMeans(psi_a(D, Dhat))
  psi <- psi(Y, D, Yhat, Dhat, thetahat)
  Bhat <- 1/n*t(psi) %*% psi
  var_thetahat <- solve(Jhat) %*% Bhat %*% t(solve(Jhat)) / n

  return(list(thetahat, var_thetahat))
}

################################################################################
#                Alternative DSR estimator - spline only                       #
################################################################################

# no covariates atm
dml_alt_nosplit_spline <- function(Y, D, S, k=100){
  require(mgcv)
  truncate <-  function(x, a, b) {
    x[x < a] <- a
    x[x > b] <- b
    x
  }
  
  psi_a <- function(D, Dhat){
    n <- nrow(D)
    p <- ncol(D)
    arr <- array(dim=c(n, p, p))
    for(i in 1:n){
      arr[i,,] <- -(D - Dhat)[i,] %*% t((D)[i,])
    }
    # -(D - Dhat) * (D-Dhat)
    return(arr)
  }
  
  psi_b <- function(Y, D, Yhat, Dhat){
    (D - Dhat) * as.numeric(Y - Yhat)
  }
  
  psi <- function(Y, D, Yhat, Dhat, thetahat){
    
    psi.a <- psi_a(D, Dhat)
    psi.b <- psi_b(Y, D, Yhat, Dhat)
    
    n <- nrow(psi.b)
    psi <- matrix(nrow=nrow(psi.b), ncol=ncol(psi.b))
    
    for(i in 1:n){
      psi[i,] <- psi.a[i,,] %*% thetahat + psi.b[i,]
    }
    return(psi)
  }
  
  n <- length(Y)
  p <- ncol(D)
  
  # estimate coefficients on full sample for speed
  # currently with no X included
  gam.dat <- data.frame(Y, D, S)
  colnames(gam.dat) <- c("Y", "D", "s1", "s2")
  fitY <- gam(Y ~ D + s(s1, s2, k=k), data=gam.dat)
  fitD <- list(length=p)
  for(j in 1:p){
    fitD[[j]] <- gam(D ~ s(s1, s2, k=k), data=gam.dat)
  }
  # Put the predictions for each fold in lists
  pred.gam.dat <- data.frame(Y, D=0, S)
  colnames(pred.gam.dat) <- c("Y", "D", "s1", "s2")
  predictionsY <- matrix(predict(fitY, newdata=pred.gam.dat), nrow=n)
  predD <- matrix(nrow=nrow(S), ncol=p)
  
  for(j in 1:p){
    predD[,j] <- predict(fitD[[j]], newdata=pred.gam.dat)
  }
  predictionsD <- predD
  
  for(j in 1:p){
    if( (min(as.matrix(D)[,j]) == 0) & (max(as.matrix(D)[,j])== 1) ){
      predictionsD[,j] <- truncate(predictionsD[,j], 0, 1)
    }
  }
  
  
  # Put the list of predictions into vectors
  Yhat <- predictionsY
  Dhat <- predictionsD
  
  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- solve(t(Vhat) %*% (D)) %*% (t(Vhat) %*% (Y - Yhat))
  # Get variance estimate
  Jhat <- colMeans(psi_a(D, Dhat))
  psi <- psi(Y, D, Yhat, Dhat, thetahat)
  Bhat <- 1/n*t(psi) %*% psi
  var_thetahat <- solve(Jhat) %*% Bhat %*% t(solve(Jhat)) / n

  return(list(thetahat, var_thetahat))
}


