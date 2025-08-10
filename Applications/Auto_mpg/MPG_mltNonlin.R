rm(list = ls())
set.seed(1234)
library(dplyr)
library(caret)
library(SeBR)

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
  
  MPG_train <- MPG[-id_save[[i]] ,]
  MPG_test <- MPG[id_save[[i]], ]
  
  n_train <- nrow(MPG_train)
  
  n_test <- nrow(MPG_test)
  
  Y_test <- MPG$MPG[id_save[[i]]]
  
  
  MPG_var <- numeric_var('MPG', support = c(-5, 5), bounds = c(-Inf, Inf))
  
  b_R <- as.basis(~cylinders + displacement + 
                    horsrpower + weight + acceleration + modelyear + origin, data = MPG_train, 
                  scale = T)
  
  
  ctmm<- ctm(response = Bernstein_basis(MPG_var, order = 4, ui = 'increasing'), 
             interacting = b_R, data = MPG_train)
  
  
  mltm <- mlt(ctmm, data = MPG_train, scale = F)
  
  pred_mltm <- predict(mltm, newdata = MPG_test[, -1], type = 'quantile', 
                       prob = c(.5, .975, .025))
  
  
  cov_mlt<- pred_mltm[-1, ]
  
  
  
  MAE_SeBR<- mean(abs(pred_mltm[1, ] - Y_test))
  
  
  CP <- as.numeric(rowSums(t(apply(cov_mlt- Y_test, 2, sign)))==0) ### coverage
  
  CI_length <- mean(abs(cov_mlt[1, ] - cov_mlt[2, ]))
  
  res <- c(MAE_SeBR, mean(CP), CI_length)
  
  outres <- rbind(outres, res)
  
  colnames(outres) <- c('MAE', 'CP', 'Length')
  
  write.table(outres, 'MPG_mlt_nonLin.txt', col.names = T, row.names = F)
  
}





