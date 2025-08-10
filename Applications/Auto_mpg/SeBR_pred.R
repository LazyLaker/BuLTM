rm(list = ls())
library(caret)
library(BuLTM)
library(SeBR)
library(dplyr)
library(doParallel)
library(parallel)

set.seed(1234)
load('MPG_aug.RData')


Csigmoid <- function(x, C=5){
  if(C <= 0){stop('C is not positive')}
  else{
    C/(1+ exp(-x))
  }
}


CIalpha <- function(x, alpha){
  if(alpha < 0.5){
    alpha = 1-alpha
  }
  return(c(quantile(x, 1-alpha), quantile(x, alpha)))
}



id_save <- list()

d0 <- NonlinAug_Data


d0$id <- 1:nrow(d0)
m <- 10
for(i in 1:m){
  train <- d0 %>% dplyr::sample_frac(0.90)
  test  <- dplyr::anti_join(d0, train, by = 'id')
  id_save[[i]] <- test$id
}

outres <- c()


for(i in 1:m){
  Y_train <- NonlinAug_Data$Y[-id_save[[i]]]
  Y_test <- NonlinAug_Data$Y[id_save[[i]]]
  
  X_train <- NonlinAug_Data[-id_save[[i]], -1]
  
  X_test <- NonlinAug_Data[id_save[[i]], -1]
  
  n_train <- nrow(X_train)
  
  n_test <- nrow(X_test)
  

  
  fit_SeBR <- sblm(Y_train, as.matrix(X_train), X_test = as.matrix(X_test))
  
  
  pred_SeBR <- fit_SeBR$fitted.values
  
  cov_SeBR<- matrix(0, n_test, 2)
  
  
  for(j in 1:n_test){
    cov_SeBR[j, ] <- CIalpha(fit_SeBR$post_ypred[, j], alpha = .05)
  }
  
  
  MAE_SeBR<- mean(abs(pred_SeBR - Y_test))
  
  
  CP <- as.numeric(rowSums(t(apply(cov_SeBR- Y_test, 1, sign)))==0) ### coverage
  
  CI_length <- mean(abs(cov_SeBR[, 1] - cov_SeBR[, 2]))
  
  res <- c(MAE_SeBR, mean(CP), CI_length)
  
  outres <- rbind(outres, res)

}



colnames(outres) <- c('MAE', 'CP', 'Length')

write.table(outres, file = 'MPG_SeBR.txt', col.names = T, row.names = F)


apply(outres[1:5, ], 2, mean)

