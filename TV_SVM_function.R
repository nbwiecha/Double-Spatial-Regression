#################################################
#   Training-validation "LS-SVM"                #
#     aka Kernel Ridge Regression               #
#     aka GP posterior mean                     #
#    with Gaussian kernel and CV'd parameters   #
#################################################

# Implementing the scheme given by Eberts and Steinwart (2013)
# Giving essentially-optimal convergence rates for the GP posterior mean
# (which is equivalent to the least-squares support vector machine (LS-SVM)
# and KRR estimators with squared error loss)

# The scheme is to use a Gaussian kernel, with cross-validated regularization
# (equiv to variance) and lengthscale parameters.

# Function will take as input the training and prediction datasets
# And the fineness of the grids used for parameter selection by CV
# And will output final predictions

# It will:
# 1. Randomly divide the training data into two splits, D1 and D2
# 2. Create a grid of (lambda, gamma) pairs (regularization, lengthscale)
# 3. For each pair, predict the values of D2 using D1 as training data
# 4. Save the MSE of the predictions on D2
# 5. Select the (lambda, gamma) pair that minimized MSE
# 6. Predict on the prediction dataset using select parameters and training data

library(fields) # for rdist() function

cv_gp_predict <- function(trainX, testX,
                          trainY, #testY,
                          n=nrow(trainX), epsilon=1/n, delta=1/n^.25,
                          silent=FALSE){
  
  # 1. Randomly split training data into D1, D2
  
  m <- floor(n/2) + 1
  
  D1.idx <- sample(1:n, m)
  D1.X <- trainX[D1.idx,]
  D1.Y <- trainY[D1.idx]
  D2.X <- trainX[-D1.idx,]
  D2.Y <- trainY[-D1.idx]

    # 2. Create a grid of (lambda, gamma) pairs
  # lambda is regularization parameter, equivalent to sigma^2 in GP regn
  # gamma is the lengthscale of the Gaussian kernel
  
  lambda.grid <- seq(0, 1, by=2*epsilon)[-1]
  gamma.grid <- seq(0, 1, by=2*delta)[-1]
 
  # 3. For each (lambda, gamma) pair, predict on D2 using D1 as training data
  d1 <- rdist(D1.X)
  d2 <- rdist(D2.X, D1.X)
  
  nrow.grid <- length(lambda.grid) * length(gamma.grid)
  
  lowest.mse <- Inf
  selected.lambda <- NA
  selected.gamma <- NA
  
  if(!silent){
    pb <- txtProgressBar(min = 0, max = nrow.grid, initial = 0, style=3) 
  }
  
  j <- 0
  
  for(gamma in gamma.grid){
    
    K1 <- exp(-d1^2/gamma^2)
    K2 <- exp(-d2^2/gamma^2)
    
    eig.K1 <- eigen(K1, symmetric=TRUE)
    d.K1 <- eig.K1$values
    Q.K1 <- eig.K1$vectors
    
    K2.Q <- K2 %*% Q.K1
    Q.y <- t(Q.K1) %*% D1.Y
    
    for(lambda in lambda.grid){
      j <- j+1
      
      if(!silent){setTxtProgressBar(pb, j)}
      # K <- K1 + lambda*m*diag(m)
      # K1.inv <- chol2inv(chol(K)) # inverse of K1
      # yhat <- K2 %*% (K1.inv %*% D1.Y)
      yhat <- K2.Q %*% ((1/(d.K1 + m*lambda))*Q.y)
      
      mse <- mean((yhat - D2.Y)^2)
      
      if(mse < lowest.mse){
        
        lowest.mse <- mse
        lambda.selected <- lambda
        gamma.selected <- gamma
        
      }
    }
    
  }
  
  if(!silent){close(pb)}
  
  # 6. Predict on test set
  d.train <- rdist(trainX)
  d.test <- rdist(testX, trainX)
  
  K.train <- exp(-d.train^2/gamma.selected^2)
  K.test <- exp(-d.test^2/gamma.selected^2)
  
  yhat.test <- K.test %*% (chol2inv(chol(K.train + n*lambda.selected*diag(n))) %*% trainY)
  return(list(yhat=yhat.test, lambda=lambda.selected, gamma=gamma.selected))
}

