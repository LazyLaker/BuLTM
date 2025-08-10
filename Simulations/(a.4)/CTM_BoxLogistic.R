rm(list = ls())
library(BuLTM)
library(dplyr)
library(distr)
library(parallel)
library(doParallel)
library(pracma)
library(SeBR)


load('BoxLogistic.RData')
load('BoxLogistic_test.RData')

colnames <- c('B1_Bayes', 'B2_Bayes', 'B3_Bayes',
              'SD1_Bayes', 'SD2_Bayes', 'SD3_Bayes',
              'CI1_Bayes', 'CI2_Bayes','CI3_Bayes', 
              'N1', 'N2', 'N3', 
              'B1_SeBR', 'B2_SeBR', 'B3_SeBR',
              'SD1_SeBRs', 'SD2_SeBR', 'SD3_SeBR',
              'CI1_SeBR', 'CI2_SeBR','CI3_SeBR', 
              'IMSEBuLTM', 'IMSESeBR', 'MAEBuLTM', 'MAESeBR')

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

###############

beta0 <- c(1, 1, 1)
beta0 <- beta0 / norm(beta0, '2')#real beta 

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
  
  real_curve <- matrix(0, n_test, length(Tgrid))
  
  for(i in 1:n_test){
    real_curve[i, ] <- 1- plogis(c(-beta0 %*% t(X_test[i, ]))+ as.numeric(sapply(grid, sign_bc)))
  }
  
  A <- BuLTM(TY, X, TMdataT$delta, probseries = seq(0, 1, .2), K = 12, iter = 1500, a = 1, eta = .01, zeta = .5)
  
  chain <-  A$chain
  MCMC_DIAG <- A$MCMC_DIAG
  hat_beta <- A$hat_beta
  sd_Bayes <- apply(chain$beta_trans, 2, sd) 
  
  chain_beta_trans <- chain$beta_trans 
  
  trans <- MCMC_DIAG %>% filter(row.names(MCMC_DIAG) 
                                %in% c("beta[1]", 
                                       "beta[2]", "beta[3]"))
  neff <- trans$n_eff
  
  
  CI <- matrix(0, p, 2)
  for (pp in 1:p) {
    CI[pp, ] <- CIalpha(chain_beta_trans[, pp], .975)
    
  }
  
  CP <- as.numeric(rowSums(t(apply(CI- beta0, 1, sign)))==0) ### coverage
  
  
  
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  
  
  pred_BuLTM <- foreach(i = 1:nrow(X_test), 
                        .combine = 'rbind', 
                        .packages = c('BuLTM')) %dopar% {predict(A, X = as.numeric(X_test[i, ]), grids = Tgrid)$estimate}
  
  stopCluster(cl)
  
  rm(cl)
  gc()
  
  
  
  
  predVMed_BuLTM <- numeric(n_test)
  
  
  for(i in 1:n_test){
    predVMed_BuLTM[i] <-  grid[which.min(abs(pred_BuLTM[i, ] - 0.5))]
  }
  
  
  
  
  fit_SeBR <- sblm(Y, as.matrix(X), X_test = as.matrix(X_test))
  
  post_SeBR <- fit_SeBR$post_ypred
  
  pred_SeBR <- matrix(0, nrow(X_test), length(grid))
  
  for(i in 1:1:nrow(X_test)){
    Fn <- ecdf(post_SeBR[, i])
    pred_SeBR[i, ] <- 1- Fn(grid)
  }
  
  predV_SeBR <- fit_SeBR$fitted.values
  
  beta_SeBR <- fit_SeBR$coefficients
  
  SeBR_sd <- apply(fit_SeBR$post_theta, 2, sd)
  
  
  CI_SeBR <- matrix(0, p, 2)
  for (pp in 1:p) {
    CI_SeBR[pp, ] <- CIalpha(fit_SeBR$post_theta, .975)
    
  }
  
  CP_SeBR <- as.numeric(rowSums(t(apply(CI_SeBR- beta0, 1, sign)))==0) ### coverage
  
  
  
  RIMSE_BuLTM <- sqrt(0.01*mean((pred_BuLTM - real_curve)^2))
  
  RIMSE_SeBR <- sqrt(0.01*mean((pred_SeBR - real_curve)^2))
  
  MAE_BuLTM <- mean(abs(predVMed_BuLTM - Y_test))
  
  MAE_SeBR <- mean(abs(predV_SeBR - Y_test))
  
  
  
  res <- c(hat_beta, sd_Bayes, CP, neff, beta_SeBR, SeBR_sd, CP_SeBR, RIMSE_BuLTM, RIMSE_SeBR, MAE_BuLTM,  MAE_SeBR)
  
  outres <- rbind(outres, res)
  
  colnames(outres) <- colnames
  
  write.table(outres, file = 'May2025_BoxLogistic.txt',
              row.names = F, col.names = T)  
  
}



