rm(list = ls())
set.seed(12300116)
library(survival)
library(xtable)
library(ggplot2)
library(MASS)
library(distr)
#######################################
inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-5)[1]
  }
}


f <- function(x){
  (0.8*x+ sqrt(x)+.825) *(0.5 *pnorm(x, 1, .3) + 0.5 * pnorm(x, 3, .3) - (0.5 *pnorm(0, 1, .3) + 0.5 * pnorm(0, 3, .3)))
}



invf <- inverse(function(x) f(x), 0, 10000000)
################################################


dat <- list()

M <- 100
n <- 200
beta0 <- c(1, 1, 1)

beta0 <- beta0 / norm(beta0, '2')#real beta 
p <- length(beta0)
rho <- .5

for(i in 1:M){
  X <- matrix(mvrnorm(n, c(0, .5), (matrix(rho, p-1, p-1)+diag(rep(1-rho, p-1)))), n)
  #const <- rt(n, 3)
  const <- rbinom(n, 1, .5)
  
  X <- cbind(const, X)
  #X_pred <- apply(X, 2, scale)
  mui <- exp(X %*% beta0)#Linear term
  #ei <- rlnorm(n, 2, .5)
  distr  <- UnivarMixingDistribution(Norm(-.5, .5), Norm(1.5, 1), mixCoeff = c(.5, .5))
  
  ei <- exp(r(distr)(n))
  
  #ei <- rf(n, 2, 1)
  #ei <- rexp(n, 2.5)
  #Ti <- sqrt(ei * mui+1/4)-1/2   ###H(t)= t(t+1)
  Ti <- as.numeric(sapply(ei*mui, invf))
  
  #Ci <- runif(n, 1, 5)  20 censoring for non AFTPO
  ##Ci <- pmin(rexp(n, 1), 5) #### 35 censoring for non AFT_PO
  
  #Ci <- pmin(rexp(n, 1), 5) #### 35 censoring for non AFT_PH
  #Ci <- pmin(rexp(n, 1), runif(n, .1, 5))  #### 55 censoring for non AFT_PH
  #Ci <- runif(n, 1, 5) #### 59Censoring for NON AFT Others
  Ci <- runif(n, 1, 3.5)
  #Ci <- runif(n, .1, 1.5)
  # Ci <- pmin(rexp(n, 1), runif(n, .1, 2))
  #Ci <- runif(n, .1, 4)
  Yi <- as.numeric(pmin(Ti, Ci)) #Observation
  #Yi <- round(Yi, 4)
  
  maxY <- max(Yi)
  delta <- as.numeric(Ti <= Ci) #Index
  TMdataT <- cbind(Yi,X, delta)
  colnames(TMdataT) <- c('Y', 'const', 'x1', 'x2', 'delta')
  TMdataT <- as.data.frame(TMdataT)
  dat[[i]] <- TMdataT
}



save(dat, file = 'Non_data_high.RData')





colMeans( do.call('rbind', dat))
