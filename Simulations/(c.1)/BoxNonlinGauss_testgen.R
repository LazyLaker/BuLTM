rm(list = ls())
set.seed(54321789)
library(survival)
library(xtable)
library(ggplot2)
library(MASS)
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

p <- 2


dat_test <- list()

for(m in 1:M){
  
  x1 <- runif(n, -2, 2)
  
  x2 <- runif(n, -2, 2)
  
  e <- rnorm(n)
  
  Z <- -x1 + pi * sin(pi*x1) + x2/2 + 15*dnorm(2*(x2-0.2)) - dnorm(x2 + 0.4) + e
  
  y <- sapply(Z, g_inv_bc)
  
  delta <- rep(1, n) #Index
  X <- cbind(x1, x2)
  TMdataT <- cbind(y,X, delta)
  colnames(TMdataT) <- c('Y', 'x1', 'x2', 'delta')
  TMdataT <- as.data.frame(TMdataT)
  write.csv(TMdataT, file = paste0('BoxNonlinGauss/test_', as.character(m), '.csv'), 
            row.names = F)
  dat_test[[m]] <- TMdataT
}

save(dat_test, file = 'BoxNonlinGauss_test.RData')


