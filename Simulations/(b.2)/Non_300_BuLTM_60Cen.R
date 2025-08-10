rm(list = ls())
library(BuLTM)
library(dplyr)
library(distr)
library(parallel)
library(doParallel)
library(SurvMetrics)

f <- function(x){
  (0.8*x+ sqrt(x)+.825) *(0.5 *pnorm(x, 1, .3) + 0.5 * pnorm(x, 3, .3) - (0.5 *pnorm(0, 1, .3) + 0.5 * pnorm(0, 3, .3)))
}


distr  <- UnivarMixingDistribution(Norm(-.5, .5), Norm(1.5, 1), mixCoeff = c(.5, .5))




load('Non_data_high.RData')
load('Non_data_high_test.RData')

colnames <- c('B1_Bayes', 'B2_Bayes', 'B3_Bayes',
              'SD1_Bayes', 'SD2_Bayes', 'SD3_Bayes',
              'CI1_Bayes', 'CI2_Bayes','CI3_Bayes', 
              'N1', 'N2', 'N3', 
              'IMSEBuLTM', 'Cindex')



beta0 <- c(1, 1, 1)
beta0 <- beta0 / norm(beta0, '2')#real beta 
grids <- seq(.1, 3.6, by = .01)

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
  
  Y_test <- TMdataT_test$Y
  X_test <- TMdataT_test[, 2:(2+p-1)]
  
  n <- length(Y)
  n_test <- length(Y_test)
  
  real_curve <- matrix(0, n_test, length(grids))
  
  for(i in 1:n_test){
    real_curve[i, ] <- 1- p(distr)(log(as.numeric(exp(-(beta0) %*% t(X_test[i, ]))) * as.numeric(sapply(grids, f))))
  }
  
  A <- BuLTMfit(Y, X, TMdataT$delta, probseries = seq(0, 1, .25), K = 12, iter = 1500, alpha = 1, eta = 1)
  
  
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
                        .packages = c('BuLTM')) %dopar% {predict(A, X = as.numeric(X_test[i, ]), grids = grids)$estimate}
  
  
  stopCluster(cl)
  
  rm(cl)
  gc()
  
  
  predVMed_BuLTM <- numeric(n_test)
  for(i in 1:n_test){
    predVMed_BuLTM[i] <-  grids[which.min(abs(pred_BuLTM[i, ] - mean(TMdataT$delta)))]
  }
  
  RIMSE_BuLTM <- sqrt(0.01*mean((pred_BuLTM - real_curve)^2))
  
  CInd_BuLTM <- Cindex(Surv(Y_test, TMdataT_test$delta), predVMed_BuLTM)
  
  res <- c(hat_beta, sd_Bayes, CP, neff, RIMSE_BuLTM, CInd_BuLTM)
  
  outres <- rbind(outres, res)
  
  colnames(outres) <- colnames
  
  
  write.table(outres, file = 'Apr182025_NonHigh_BuLTM.txt',
              row.names = F, col.names = T)  
}


