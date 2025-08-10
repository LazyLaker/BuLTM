rm(list = ls())
set.seed(24250102)
library(survival)
library(xtable)
library(ggplot2)
library(MASS)

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


M <- 100
n <- 20
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

dat_test <- list()

for(m in 1:M){
  X <- mvrnorm(n, c(0, 0, 0), Sigma = Sigma)
  
  Z <- X %*% beta0 + log(rexp(n))
  
  y <- sapply(Z, g_inv_bc)
  
  delta <- rep(1, n) #Index
  TMdataT <- cbind(y,X, delta)
  colnames(TMdataT) <- c('Y', 'x1', 'x2', 'x3', 'delta')
  TMdataT <- as.data.frame(TMdataT)
  write.csv(TMdataT, file = paste0('BoxPH/test_', as.character(m), '.csv'), 
            row.names = F)
  
  dat_test[[m]] <- TMdataT
}

save(dat_test, file = 'BoxPH_test.RData')


