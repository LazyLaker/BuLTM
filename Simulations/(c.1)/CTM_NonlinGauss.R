rm(list = ls())
library(BuLTM)
library(dplyr)
library(parallel)
library(doParallel)
library(pracma)
library(SeBR)
library(splines2)

load('BoxNonlinGauss.RData')
load('BoxNonlinGauss_test.RData')

colnames <- c('IMSEBuLTM', 'IMSESeBR', 'MAEBuLTM', 'MAESeBR')

outres <- c()

####### functions

Csigmoid <- function(x, C=5){
  if(C <= 0){stop('C is not positive')}
  else{
    C/(1+ exp(-x))
  }
}

Clogit <- function(x, C= 5){
  log((x/C)/(1-x/C))
}


sign_bc <- function(s, lambda = 2){
  if(lambda == 0) {
    # Inverse log-transformation
    exp(s)
  } else{
    (sign(s)* abs(s)^(lambda) - 1)/lambda
  }
}

CIalpha <- function(x, alpha){
  if(alpha < 0.5){
    alpha = 1-alpha
  }
  return(c(quantile(x, 1-alpha), quantile(x, alpha)))
}


outres <- c()


for(j in 1:100){
  
  TMdataT <- dat[[j]]
  
  TMdataT_test <- dat_test[[j]]
  
  
  p <- ncol(TMdataT) - 2
  
  
  Y <- TMdataT$Y
  X <- TMdataT[, 2:(2+p-1)]
  
  TY <- sapply(Y, Csigmoid)
  
  Y_test <- TMdataT_test$Y
  X_test <- TMdataT_test[, 2:(2+p-1)]
  
  n <- length(Y)
  n_test <- length(Y_test)
  
  grid <- seq(min(Y)-0.1, max(Y)+0.1, by = .01)
  
  Tgrid <- sapply(grid, Csigmoid)
  
  S1 <- bSpline(X$x1, knots = quantile(X$x1, seq(0.1, .9, .1)), Boundary.knots = c(-2, 2))
  
  S2 <- bSpline(X$x2, knots = quantile(X$x2, seq(0.1, .9, .1)), Boundary.knots = c(-2, 2))
  
  TX <- cbind(S1, S2)
  
  delta <- TMdataT$delta
  
  model <- BuLTM(TY, TX, delta = delta, probseries = seq(0, 1, .2), K = 12, a = 1, eta = .01, zeta = .5)
  
  
  TX_test_1 <- predict(S1, X_test[, 1])
  TX_test_2 <- predict(S2, X_test[, 2])
  TX_test <- cbind(TX_test_1, TX_test_2)
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  pred_BuLTM <- foreach(i = 1:nrow(X_test), 
                        .combine = 'rbind', 
                        .packages = c('BuLTM')) %dopar% {
                          predict(model, X = as.numeric(TX_test[i, ]), grids = Tgrid)$estimate}
  
  stopCluster(cl)
  
  rm(cl)
  gc()
  
  
  real_curve <- matrix(0, n_test, length(Tgrid))
  
  for(i in 1:n_test){
    x1i <- X_test[i, 1]
    x2i <- X_test[i, 2]
    zi <- -x1i + pi * sin(pi*x1i) + x2i/2 + 15*dnorm(2*(x2i-0.2)) - dnorm(x2i + 0.4)
    real_curve[i, ] <- 1- pnorm(-zi+ as.numeric(sapply(grid, sign_bc)))
  }
  
  
  predVMed_BuLTM <- numeric(n_test)
  
  
  for(i in 1:n_test){
    predVMed_BuLTM[i] <-  grid[which.min(abs(pred_BuLTM[i, ] - 0.5))]
  }
  
  
  
  RIMSE_BuLTM <- sqrt(0.01*mean((pred_BuLTM - real_curve)^2))
  
  MAE_BuLTM <- mean(abs(predVMed_BuLTM - Y_test))
  
  
  
  
  fit_SeBR <- sblm(Y, as.matrix(TX), X_test = as.matrix(TX_test))
  
  post_SeBR <- fit_SeBR$post_ypred
  
  pred_SeBR <- matrix(0, nrow(X_test), length(grid))
  
  for(i in 1:1:nrow(X_test)){
    Fn <- ecdf(post_SeBR[, i])
    pred_SeBR[i, ] <- 1- Fn(grid)
  }
  
  predV_SeBR <- fit_SeBR$fitted.values
  
  RIMSE_SeBR <- sqrt(0.01*mean((pred_SeBR - real_curve)^2))
  
  MAE_SeBR <- mean(abs(predV_SeBR - Y_test))
  
  
  res <- c(RIMSE_BuLTM, RIMSE_SeBR, MAE_BuLTM, MAE_SeBR)
  
  outres <- rbind(outres, res)
  
  colnames(outres) <- colnames
  
  write.table(outres, file = 'May2025_NonlinGauss.txt',
              row.names = F, col.names = T)  
  
  
  
}


