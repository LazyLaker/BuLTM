rm(list = ls())
set.seed(2020666)
library(survival)
library(xtable)
library(ggplot2)
library(MASS)
library(distr)

library(SeBR)
lambda <- .5

g_inv_bc = function(s, lambda = 2) {
  if(lambda == 0) {
    # Inverse log-transformation
    exp(s)
  } else {
    # Inverse (signed) Box-Cox-transformation
    sign(lambda*s + 1)*abs(lambda*s+1)^(1/lambda)
  }
}


M <- 300
n <- 200
beta0 <- c(1, 1, 1)

beta0 <- beta0 / norm(beta0, '2')#real beta 
p <- length(beta0)
rho <- .75


Sigma <- matrix(0, p, p)



for(i in 1:p){
  for(j in 1:p){
    Sigma[i, j] <- .75^(abs(i-j))
  }
}

dat <- list()

distr  <- UnivarMixingDistribution(Norm(-.5, .5), Norm(1.5, 1), mixCoeff = c(.5, .5))

for(m in 1:M){
  X <- mvrnorm(n, c(0, 0, 0), Sigma = Sigma)
  
  e <- r(distr)(n)
  
  Z <- X %*% beta0 + e
  
  y <- sapply(Z, g_inv_bc)
  
  delta <- rep(1, n) #Index
  TMdataT <- cbind(y,X, delta)
  colnames(TMdataT) <- c('Y', 'x1', 'x2', 'x3', 'delta')
  TMdataT <- as.data.frame(TMdataT)
  dat[[m]] <- TMdataT
}

save(dat, file = 'BoxMix.RData')

