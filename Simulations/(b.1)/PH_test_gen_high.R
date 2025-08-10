rm(list = ls())

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
  (0.8*x+ sqrt(x)+.825) *(0.5 *pnorm(x, .5, .2) + 0.5 * pnorm(x, 2.5, .3) - (0.5 *pnorm(0, 0.5, .2) + 0.5 * pnorm(0, 2.5, .3)))
}


invf <- inverse(function(x) f(x), 0, 10000000)
################################################

set.seed(2025666)

dat_test <- list()

M <- 100
n <- 20
beta0 <- c(1, 1, 1)

beta0 <- beta0 / norm(beta0, '2')#real beta 
p <- length(beta0)
rho <- .5

for(i in 1:M){
  X <- matrix(mvrnorm(n, c(0, 0), (matrix(rho, p-1, p-1)+diag(rep(1-rho, p-1)))), n)
  #const <- rt(n, 3)
  const <- rbinom(n, 1, .5)
  
  X <- cbind(const, X)
  #X_pred <- apply(X, 2, scale)
  mui <- exp(X %*% beta0)#Linear term
  #ei <- rlnorm(n, 2, .5)
  
  ei <- rexp(n, 1)
  # distr  <- UnivarMixingDistribution(Exp(2.5),  Exp(1/3), mixCoeff = c(.55, .45))
  # ei <- r(distr)(n)
  #ei <- exp(rlogis(n))
  
  #ei <- rf(n, 2, 1)
  #ei <- rexp(n, 2.5)
  #Ti <- sqrt(ei * mui+1/4)-1/2   ###H(t)= t(t+1)
  Ti <- as.numeric(sapply(ei*mui, invf))
  
  #Ci <- runif(n, 1, 5)  20 censoring for non AFTPO
  ##Ci <- pmin(rexp(n, 1), 5) #### 35 censoring for non AFT_PO
  
  #Ci <- pmin(rexp(n, 1), 5) #### 35 censoring for non AFT_PH
  #Ci <- pmin(rexp(n, 1), runif(n, .1, 5))  #### 55 censoring for non AFT_PH
  #Ci <- runif(n, 1, 5) #### 59Censoring for NON AFT Others
  Ci <- pmin(rexp(n, 1), 2.5)
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
  dat_test[[i]] <- TMdataT
}



save(dat_test, file = 'PH_data_high_test.RData')

colMeans( do.call('rbind', dat_test))

