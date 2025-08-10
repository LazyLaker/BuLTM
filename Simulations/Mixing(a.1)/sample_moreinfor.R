rm(list = ls())
library(BuLTM)
library(bayesplot)
set.seed(1234)

load('BoxGauss.RData')


Csigmoid <- function(x, C=5){
  if(C <= 0){stop('C is not positive')}
  else{
    C/(1+ exp(-x))
  }
}


TMdataT <- dat[[1]]

p <- ncol(TMdataT) - 2


Y <- TMdataT$Y
X <- TMdataT[, 2:(2+p-1)]

TY <- sapply(Y, Csigmoid)

n <- length(Y)

delta <- TMdataT$delta



eta = .01
zeta = .5

model <- BuLTM(TY, X, delta = delta, K = 12, eta = eta, iter = 1000, 
               zeta = zeta, probseries = seq(0, 1, .2), max_treedepth = 10)


SuffitInf(model)

mcmc_trace(model$stan, pars = 'lp__')

A <- model$stan


mcmc_trace(A, pars = c('beta[1]', 'beta[2]', 'beta[3]')) 

### (384, 1.003), (480, 1.006), (416. 1.002)



