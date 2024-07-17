# Generic DML Function for Partially Linear Model
# For Two-Stage Estimators for Spatial Confounding
# Code by Nate Wiecha, North Carolina State University

# predict_fn will need to be specified, it makes predictions for Y and D
# Y is response, D treatment (vector), S spatial locations, X covariate
# K is number of folds, folds is vector of fold assignments

# predict_fn will take as input trainY, trainD, trainS, testS, (trainX, testX)
# and output a list with components predY, predD

dml_lhat <- function(Y, D, S, K, predict_fn, ...){
  
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
  folds <- sample(1:K, size=n, replace=TRUE)
  
  # Put the predictions for each fold in lists
  predictionsY <- predictionsD <- list(length=K)
  for(k in 1:K){
    trainS <- S[folds!=k,]
    trainD <- D[folds!=k,]
    trainY <- Y[folds!=k]
    testS <- S[folds==k,]
    
    predictions <- predict_fn(trainS, trainY, trainD, testS, ...)
    predictionsY[[k]] <- predictions$predY
    predictionsD[[k]] <- predictions$predD
    for(j in 1:p){
      if( (min(as.matrix(trainD)[,j]) == 0) & (max(as.matrix(trainD)[,j])== 1) ){
        predictionsD[[k]][,j] <- truncate(predictionsD[[k]][,j], 0, 1)
      }
    }
  }
  
  # Put the list of predictions into vectors
  Yhat <- rep(NA, length(Y))
  Dhat <- matrix(nrow=nrow(D), ncol=p)
  for(k in 1:K){
    Yhat[folds==k] <- predictionsY[[k]]
    Dhat[folds==k,] <- predictionsD[[k]]
  }
  
  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- lm((Y-Yhat) ~ (Vhat)-1)$coefficients
  # Get variance estimate
  Jhat <- colMeans(psi_a(D, Dhat))
  psi <- psi(Y, D, Yhat, Dhat, thetahat)
  Bhat <- 1/n*t(psi) %*% psi
  var_thetahat <- solve(Jhat) %*% Bhat %*% t(solve(Jhat)) / n

  return(list(thetahat, var_thetahat))
}

################################################################################
#                Another function that just includes covariates in             #
#                        predictive model                                      #
################################################################################

dml_lhat_covariates <- function(Y, D, S, X, K, predict_fn, ...){
  
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
  folds <- sample(1:K, size=n, replace=TRUE)
  
  # Put the predictions for each fold in lists
  predictionsY <- predictionsD <- list(length=K)
  for(k in 1:K){
    trainS <- S[folds!=k,]
    trainD <- D[folds!=k,]
    trainX <- X[folds!=k,]
    trainY <- Y[folds!=k]
    testS <- S[folds==k,]
    testX <- X[folds==k,]
    
    predictions <- predict_fn(trainS=trainS, trainY=trainY, trainD=trainD, 
                              testS=testS, trainX=trainX, testX=testX, ...)
    predictionsY[[k]] <- predictions$predY
    predictionsD[[k]] <- predictions$predD
    for(j in 1:p){
      if( (min(as.matrix(trainD)[,j]) == 0) & (max(as.matrix(trainD)[,j])== 1) ){
        predictionsD[[k]][,j] <- truncate(predictionsD[[k]][,j], 0, 1)
      }
    }
  }
  
  # Put the list of predictions into vectors
  Yhat <- rep(NA, length(Y))
  Dhat <- matrix(nrow=nrow(D), ncol=p)
  for(k in 1:K){
    Yhat[folds==k] <- predictionsY[[k]]
    Dhat[folds==k,] <- predictionsD[[k]]
  }
  
  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- lm((Y-Yhat) ~ (Vhat)-1)$coefficients
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

dml_alt <- function(Y, D, S, X, K){
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
  folds <- sample(1:K, size=n, replace=TRUE)
  
  # estimate coefficients on full sample for speed
  fitY <- fit_model(y=Y, locs=S, X=cbind(1, D, X),silent=TRUE)
  fitD <- list(length=p)
  for(j in 1:p){
    fitD[[j]] <- fit_model(y=D[,j], locs=S, X=cbind(rep(1,n),X), silent=TRUE)
  }
  # Put the predictions for each fold in lists
  predictionsY <- predictionsD <- list(length=K)
  for(k in 1:K){
    trainS <- S[folds!=k,]
    trainD <- D[folds!=k,]
    trainX <- X[folds!=k,]
    trainY <- Y[folds!=k]
    testS <- S[folds==k,]
    testX <- X[folds==k,]

    predictionsY[[k]] <- predictions(fitY, locs_pred=testS, 
                                     X_pred=cbind(1, matrix(0, nrow=nrow(testS), ncol=p), testX))
    
    predD.k <- matrix(nrow=nrow(testS), ncol=p)
    for(j in 1:p){
      
      predD.k[,j] <- predictions(fitD[[j]], locs_pred=testS,
                                 X_pred=cbind(rep(1, nrow(testS)), testX))
    }
    
    predictionsD[[k]] <- predD.k
    for(j in 1:p){
      if( (min(as.matrix(trainD)[,j]) == 0) & (max(as.matrix(trainD)[,j])== 1) ){
        predictionsD[[k]][,j] <- truncate(predictionsD[[k]][,j], 0, 1)
      }
    }
  }
  
  # Put the list of predictions into vectors
  Yhat <- rep(NA, length(Y))
  Dhat <- matrix(nrow=nrow(D), ncol=p)
  for(k in 1:K){
    Yhat[folds==k] <- predictionsY[[k]]
    Dhat[folds==k,] <- predictionsD[[k]]
  }
  
  # Get theta_hat
  Vhat <- D - Dhat
  
  thetahat <- solve(t(Vhat) %*% (D)) %*% (t(Vhat) %*% (Y - Yhat))
  # Get variance estimate
  Jhat <- colMeans(psi_a(D, Dhat))
  psi <- psi(Y, D, Yhat, Dhat, thetahat)
  Bhat <- 1/n*t(psi) %*% psi
  var_thetahat <- solve(Jhat) %*% Bhat %*% t(solve(Jhat)) / n
  # se_thetahat <- sqrt(sigma2_hat/n)
  
  return(list(thetahat, var_thetahat))
}
