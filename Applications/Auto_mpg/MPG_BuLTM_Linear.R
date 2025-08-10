rm(list = ls())
set.seed(1234)
library(dplyr)
library(caret)
library(BuLTM)

load('MPG.RData')
Csigmoid <- function(x, C=5){
  if(C <= 0){stop('C is not positive')}
  else{
    C/(1+ exp(-x))
  }
}


id_save <- list()

d0 <- MPG


d0$id <- 1:nrow(d0)
m <- 10
for(i in 1:m){
  train <- d0 %>% dplyr::sample_frac(0.90)
  test  <- dplyr::anti_join(d0, train, by = 'id')
  id_save[[i]] <- test$id
}

outres <- c()


for(i in 1:m){
  Y_train <- MPG$MPG[-id_save[[i]]]
  Y_test <- MPG$MPG[id_save[[i]]]
  
  X_train <- MPG[-id_save[[i]], -1]
  
  X_test <- MPG[id_save[[i]], -1]
  
  n_train <- nrow(X_train)
  
  n_test <- nrow(X_test)
  TY <- sapply(Y_train, FUN = function(x) Csigmoid(x))
  
  
  A <- BuLTMfit(TY, X_train, delta = rep(1, n_train), probseries = seq(0, 1, .2), K = 12, iter = 1500, alpha = 1, eta = 1)
  
  grid <- seq(min(Y_train)-.5, max(Y_train)+.5, by = .01)
  
  Tgrid <- sapply(grid, function(x) Csigmoid(x))
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  
  
  pred_BuLTM <- foreach(i = 1:n_test, 
                        .combine = 'rbind', 
                        .packages = c('BuLTM')) %dopar% {predict(A, X = as.numeric(X_test[i, ]), grids = Tgrid)$estimate}
  
  
  stopCluster(cl)
  
  
  rm(cl)
  gc()
  
  
  predVMed_BuLTM <- numeric(n_test)
  
  cov_BuLTM <- matrix(0, n_test, 2)
  
  
  for(i in 1:n_test){
    predVMed_BuLTM[i] <-  grid[which.min(abs(pred_BuLTM[i, ] - 0.5))]
    cov_BuLTM[i, 1] <- grid[which.min(abs(pred_BuLTM[i, ] - 0.025))]
    cov_BuLTM[i, 2] <- grid[which.min(abs(pred_BuLTM[i, ] - 0.975))]
  }
  
  
  MAE_BuLTM <- mean(abs(predVMed_BuLTM - Y_test))
  
  CP <- as.numeric(rowSums(t(apply(cov_BuLTM- Y_test, 1, sign)))==0) ### coverage
  
  CI_length <- mean(abs(cov_BuLTM[, 1] - cov_BuLTM[, 2]))
  
  res <- c(MAE_BuLTM, mean(CP), CI_length)
  
  outres <- rbind(outres, res)
  
  colnames(outres) <- c('MAE', 'CP', 'Length')
  
  write.table(outres, 'MPG_BuLTM_Lin.txt', col.names = T, row.names = F)
  
}

