rm(list = ls())
set.seed(1234)
library(dplyr)
library(caret)

library(extRemes)


rm(list = ls())
packages <- c("stringr","MASS","doParallel", "loo",  "BayesX","microbenchmark", "Rcpp", "RcppArmadillo", "splines", "mgcv", "Matrix", "MCMCpack", "sdPrior", "R2BayesX",
              "RhpcBLASctl", "scam", "bamlss", "rlang", "scoringutils", "scales")



source("code/helpers.R")

load_inst_packages(packages)

# restrict parallel computing 
# library("RhpcBLASctl")
# omp_set_num_threads(1)
# blas_set_num_threads(1)

# code provided by hothorn et al. (2017)
source("code/mlt_sims/simfuns.R")

# posteriors and gradients in rcpp
sourceCpp("code/rcpp/posterior_grad_xx2.cpp")



source("code/bctm_utils.R")
source("code/bctm_design_funs2.R")
source("code/bctm_design.R")
source("code/bctm_fun.R")

source("code/nuts/nuts_utils.R")
source("code/nuts/nuts_lin.R")
source("code/nuts/adnuts_helper.R")


load('MPG.RData')

library(mlt)

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
  
  
  data_train <- cbind(Y_train, X_train)
  
  colnames(data_train)[1] <- 'y'
  
  grid <- seq(min(Y_train)-3, max(Y_train)+3, by = .01)
  
  varn <- colnames(X_train)
  
  
  
  
  f <- paste("y", "~", "hy_lin(y) +", paste0("hx_lin(" , varn,")", 
                                             collapse=" + "))

  mod <- paste("~", "y+", paste0(varn,collapse=" + ")) %>% as.formula
  
  
  res <- bctm(      
    as.formula(f),
    family = "gauss", data= data_train,
    iterations = 4000, warmup = 2000, burnin = 2000,
    hyperparams=list(a=1, b=.001), nuts_settings=list(adapt_delta = 0.8, max_treedepth=16), seed = 1234)
  
  
  exp_ident <-res$model$exp_ident
  
  # beta samples
  betas <- res$samples$beta[(2000+1):4000, ]
  
  #beta_tilde samples
  bts <- betas
  bts[,exp_ident] <- exp(bts[,exp_ident])
  
  bt <- colMeans(bts)
  
  
  
  
  BCTM_pred <- matrix(0, nrow(X_test), length(grid))
  
  
  data_test <- cbind(Y_test, X_test)
  
  colnames(data_test)[1] <- 'y'
  
  
  for(i in 1:n_test){
    for(j in 1:length(grid)){
      xsmatb<- model.matrix(mod,  data=data.frame(y = grid[j], X_test[i, varn])) %*%t(bts)
      BCTM_pred[i,j] <- 1-apply(pnorm(xsmatb), 1, median)
    }
  }
  
  
  predVMed_BCTM <- numeric(n_test)
  
  CI_BCTM <- matrix(0, n_test, 2)
  
  for(i in 1:n_test){
    predVMed_BCTM[i] <-  grid[which.min(abs(BCTM_pred[i, ] - 0.5))]
    CI_BCTM[i, 1] <- grid[which.min(abs(BCTM_pred[i, ] - .025))]
    CI_BCTM[i, 2] <- grid[which.min(abs(BCTM_pred[i, ] - .975))]
    
  }
  
  MAE_BCTM <- mean(abs(Y_test - predVMed_BCTM))
  
  
  CP <- as.numeric(rowSums(t(apply(CI_BCTM- Y_test, 1, sign)))==0) ### coverage
  
  CI_length <- mean(abs(CI_BCTM[1, ] - CI_BCTM[2, ]))
  
  res <- c(MAE_BCTM, mean(CP), CI_length)
  
  outres <- rbind(outres, res)
  
  colnames(outres) <- c('MAE', 'CP', 'Length')
  
  write.table(outres, 'MPG_BCTM_Lin.txt', col.names = T, row.names = F)
  
}





