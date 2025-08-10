rm(list = ls())
library(coda)
library(survival)
library(spBayesSurv)
library(doParallel)
library(fields)
library(rstan)
library(dplyr)
library(splines2)
library(xtable)
library(AUC)
library(mlt)
library(survivalROC)
library(BuLTM)
library(caret)
library(dplyr)
library(SurvMetrics)
library(TransModel)
library(foreach)

CIalpha <- function(x, alpha){
  if(alpha < 0.5){
    alpha = 1-alpha
  }
  return(c(quantile(x, 1-alpha), quantile(x, alpha)))
}

set.seed(12300116)
############## Read table ##############

d <- read.csv('heart_failure.csv', header = T)
n <- nrow(d)
d$time <- d$time / 30 ### convert days to months
d$age <- d$age / 100 #### age over 10

d$platelets <- log(d$platelets) / 10 #### log of platelets/mL / 10
d$serum_sodium <- d$serum_sodium / 1000 ##### Eq/L
d$ejection_fraction <- scale(d$ejection_fraction)  ###### proportion of blood leaving the heart
d$creatinine_phosphokinase <- scale(log(d$creatinine_phosphokinase)) ### log phosphokinase

# 


colnames <- colnames(d)
Y <- d$time
X <- d[, -which(colnames == 'time')]
delta <- d$DEATH_EVENT
X <- X[, -ncol(X)]
colnamesX <- colnames(X)
X <- as.matrix(X)
p <- ncol(X)




S <- Surv(d$time, event = d$DEATH_EVENT, type = "right")
nburn=4000; nsave=6000; nskip=20; niter = nburn+nsave;
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(maxL=15, a0=1, b0=1, M=10, q=.9);

# dataframX <- data.frame(X)

# resPO <- survregbayes(S ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11, data = dataframX, survmodel = "PO", prior = prior, mcmc = mcmc,
#                       dist = "loglogistic")
# beta_PO <- -resPO$coefficients[1:ncol(X)]
# 
# resPH <-  survregbayes(S ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11, data = dataframX, survmodel = "PH", prior = prior, mcmc = mcmc,
#                        dist = "loglogistic")
# 
# beta_PH <- -resPH$coefficients[1:ncol(X)]
# 
# # 
# resAFT <- survregbayes(S ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11, data = dataframX, survmodel = "AFT", prior = prior, mcmc = mcmc,
#                        dist = "loglogistic")
# beta_AFT <- -resAFT$coefficients[1:ncol(X)]
# # 
# # 
# # norm_PH <- norm(beta_PH, '2')
# # # norm_PO <- norm(beta_PO, '2')
# # # norm_AFT <- norm(beta_AFT, '2')
# # # hat_beta_trans <- hat_beta * norm_PH / norm(hat_beta, '2')
# # 
#  which.min(c(resPO$DIC, resPH$DIC, resAFT$DIC))
# # 
#  which.max(c(sum(log(resPO$cpo)), sum(log(resPH$cpo)), sum(log(resAFT$cpo))))
# 
# beta_trans <- hat_beta / norm(hat_beta, '2') * norm(beta_PH, '2')
# 
# output <- data.frame(hat_beta, beta_PH)
# names(output) <- c('Proposed', 'PH')
# rownames(output) <- colnames(d[, 1:11])
# 
# 
# newdata <- apply(X, 2, mean)
# grids <- seq(0, max(t), .1)
# gridSpline <- t(iSpline(grids, intercept = F, df = 3, knots = knots, Boundary.knots = c(0, max(grids))))
# 
# 
# plot(grids, npBtrans, type = 'l', ylim = c(0, 1))
# newdata <- data.frame(t(newdata))
# names(newdata) <- colnames(dataframX)
# spPH <- plot.survregbayes(resPH, tgrid = grids, xnewdata = newdata, PLOT = F)
# spPO <- plot.survregbayes(resPO, tgrid = grids, xnewdata = newdata, PLOT = F)
# spAFT <- plot.survregbayes(resAFT, tgrid = grids, xnewdata = newdata, PLOT = F)
# 
# lines(grids, spPH$Shat, lty = 2, col = 'red')
#  lines(grids, spPO$Shat, lty = 4, col = 'blue')
#  lines(grids, spAFT$Shat, lty = 5, col = 'yellow')
# 
# 
# 
# # 
# S_AUC <- matrix(0, 4, 9)
# 
# for (i in 1:9) {
#   S_AUC[1, i] <- survivalROC(TMdataT$time, delta, marker = -X %*% hatbeta, predict.time = i, method = 'NNE', span=0.25*n^(-0.20))$AUC
#   S_AUC[2, i] <- survivalROC(TMdataT$time, delta, marker = -X %*% beta_PH, predict.time = i, method = 'NNE', span=0.25*n^(-0.20))$AUC
#   S_AUC[3, i] <- survivalROC(TMdataT$time, delta, marker = -X %*% beta_PO, predict.time = i, method = 'NNE', span=0.25*n^(-0.20))$AUC
#   S_AUC[4, i] <- survivalROC(TMdataT$time, delta, marker = -X %*% beta_AFT, predict.time = i, method = 'NNE', span=0.25*n^(-0.20))$AUC
# 
# }
# 
# plot(1:9, S_AUC[1, ], type = 'l', col = 'red', xlab = 'Follow-up time (months)', ylab = 'AUC', ylim = c(0, 1))
# lines(1:9, S_AUC[2, ], type = 'l', col = 'blue', lty = 2)
# lines(1:9, S_AUC[3, ], type = 'l', col = 'green', lty = 3)
# lines(1:9, S_AUC[4, ], type = 'l', col = 'black', lty = 4)
# 
# 
# 

### splite data
d2 <- d
d2$id <- 1:nrow(d2)


id_save <- list()

m <- 10
for(i in 1:m){
  train <- d2 %>% dplyr::sample_frac(0.90)
  test  <- dplyr::anti_join(d2, train, by = 'id')
  id_save[[i]] <- test$id
}

### prediction evaluation

outres <- c()

res_name <- c('CindexBuLTM', 'IBSBuLTM', 
              'CindexPH', 'IBSPH', 
              'CindexTM', 'IBSTM', 
              'Cindmlt', 'IBSmlt')

# CindexBuLTM <- numeric(10)
# CindexPH <- numeric(10)
# CindexTM <- numeric(10)
# Cindmlt <- numeric(10)
# 
# IBSBuLTM <- numeric(10)
# IBSPH <- numeric(10)
# IBSTM <- numeric(10)
# IBSmlt <- numeric(10)


maxtime <- max(d$time)
grids <- seq(.1, max(d$time)+.05, by = .02)

len <- length(grids)


for(j in 1:10){
  
  test <- d[id_save[[j]], ]
  train <- d[-id_save[[j]], ]
  
  trainX <- train[, 1:p]
  testX <- test[, 1:p]
  
  colnames <- c()
  for(cc in 1:ncol(trainX)){
    colnames <- append(colnames, paste0('x', as.character(cc)))
  }
  names(trainX) <- colnames
  names(testX) <- colnames
  
  
  fit <- BuLTM(Y = train$time, X= trainX, delta = train$DEATH_EVENT, 
               K = 12, a = 1, probseries = seq(0, 1, by = .2), 
               zeta = .5, eta = .01, max_treedepth = 10)
  

  S_train <- Surv(train$time, train$DEATH_EVENT, type = "right")
  
  fitPH <- survregbayes(S_train ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11, data = trainX, survmodel = "PH", prior = prior, mcmc = mcmc,
                        dist = "loglogistic")
  
  fitTM <- TransModel(S_train ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11, data = trainX, r=0)
  
  
  trainX2 <- trainX
  
  Var_y <- numeric_var("y", support = c(0, maxtime), bounds = c(0, Inf))
  
  trainX2$y <- with(trainX2, Surv(train$time, train$DEATH_EVENT))
  
  trainX2$x3 <- as.numeric(trainX2$x3)
  
  trainX2$x5 <- as.numeric(trainX2$x5)
  
  B_y <- Bernstein_basis(var = Var_y, order = 10, ui = "increasing")
  
  fm <- trainX2$y~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11
  
  ctmm <- ctm(B_y, shifting = fm[-2L], 
              data = trainX2,todistr = "MinExtrVal")
  
  res_mlt <- mlt(ctmm, data = trainX2)
  
  
  
  BuLTM_med <- numeric(nrow(test))
  PH_med <- numeric(nrow(test))
  TM_med <- numeric(nrow(test))
  mlt_med <- numeric(nrow(test))
  
  
  sp_mat_BuLTM <- matrix(0, nrow(test), len)
  sp_mat_PH <- matrix(0, nrow(test), len)
  TM_mat_PH <- matrix(0, nrow(test), len)
  
  mlt_mat_PH <-  t(predict(res_mlt, newdata = testX, type = 'survivor', K = len))
  mlt_med <- apply(mlt_mat_PH, 1, FUN = function(x) grids[which.min(abs(x - 0.7))])
  
  
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  
  
  sp_mat_BuLTM <- foreach(i = 1:nrow(testX), 
                        .combine = 'rbind', 
                        .packages = c('BuLTM')) %dopar% {predict(fit, X = as.numeric(testX[i, ]), grids = grids)$estimate}
  
  
  stopCluster(cl)
  
  rm(cl)
  gc()
  
  BuLTM_med <- apply(sp_mat_BuLTM, 1, FUN = function(x) grids[which.min(abs(x - 0.7))])
  
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  sp_mat_PH <- foreach(i = 1:nrow(testX), 
                          .combine = 'cbind', 
                          .packages = c('spBayesSurv')) %dopar% {
                            plot.survregbayes(fitPH, tgrid = grids, xnewdata = testX[i, 1:p], PLOT = F)$Shat
                          }
  # 
  TM_mat_PH <- foreach(i = 1:nrow(testX),
                       .combine = 'rbind',
                       .packages = c('TransModel')) %dopar% {
                         predict(fitTM, new.time = grids, newdata = testX[i, 1:p])$survival
                       }

  stopCluster(cl)
  
  rm(cl)
  gc()
  
  sp_mat_PH <- t(sp_mat_PH)
  
  PH_med <- apply(sp_mat_PH, 1, FUN = function(x) grids[which.min(abs(x - 0.7))])
  
  TM_med <- apply(TM_mat_PH, 1, FUN = function(x) grids[which.min(abs(x - 0.7))])
  
  
  CindexBuLTM <- Cindex(Surv(test$time, test$DEATH_EVENT), BuLTM_med)
  CindexPH <- Cindex(Surv(test$time, test$DEATH_EVENT), PH_med)
  CindexTM <- Cindex(Surv(test$time, test$DEATH_EVENT), TM_med)
  Cindmlt <- Cindex(Surv(test$time, test$DEATH_EVENT), mlt_med)
  
  
  # BrierBuLTM[j] <- Brier(Surv(test$time, test$DEATH_EVENT), BuLTM_med)
  # BrierPH[j] <- Brier(Surv(test$time, test$DEATH_EVENT), PH_med)
  
  
  # 
  IBSBuLTM <- IBS(Surv(test$time, test$DEATH_EVENT), sp_mat_BuLTM, grids)
  IBSPH <- IBS(Surv(test$time, test$DEATH_EVENT), sp_mat_PH, grids)
  IBSTM <- IBS(Surv(test$time, test$DEATH_EVENT), TM_mat_PH, grids)
  IBSmlt <- IBS(Surv(test$time, test$DEATH_EVENT), mlt_mat_PH, grids)
  
  
  res <- c(CindexBuLTM, IBSBuLTM, CindexPH, IBSPH, 
           CindexTM, IBSTM,  Cindmlt, IBSmlt)
  
  outres <- rbind(outres, res)
  
}

colnames(outres) <- res_name

# 
# data_BS = data.frame('IBS' = c(IBSBuLTM, IBSPH, IBSTM),
#                      'model' = c(rep('BuLTM', 10), rep('spBayesSurv', 10), rep('TransModel', 10)))
# 
# ggplot(data_BS, aes(x = model, y = IBS, fill = model)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#FFBBCC", "#88CCFF", '#99FFCC')) +
#   theme(legend.position="none") +
#   theme(axis.title.x = element_blank())
# 
# 
# data_CI = data.frame('Cindex' = c(CindexBuLTM, CindexPH, CindexTM),
#                      'model' = c(rep('BuLTM', 10), rep('spBayesSurv', 10), rep('TransModel', 10))) 
# 
# ggplot(data_CI, aes(x = model, y = Cindex, fill = model)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#FFBBCC", "#88CCFF", '#99FFCC')) +
#   theme(legend.position="none") +
#   theme(axis.title.x = element_blank())
# 
# 
# data_CI <- read.table('CindexAll.txt', header = T)
# 
# data_BS <- read.table('IBSAll.txt', header = T) 
