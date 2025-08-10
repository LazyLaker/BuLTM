rm(list = ls())
library(survival)
library(spBayesSurv)
library(BuLTM)
library(xtable)
library(AUC)
library(HDInterval)
library(loo)
library(survivalROC)
library(bayesplot)
library(ggplot2)
library(latex2exp)
library(survAUC)
library(latex2exp)
library(dplyr)
library(SurvMetrics)
library(caret)
library(TransModel)


########################
#data("veteran")
d <- veteran


set.seed(19671230)
d <- veteran
d2 <- d
d2$id <- 1:nrow(d2)


id_save <- list()

m <- 10
for(i in 1:m){
  train <- d2 %>% dplyr::sample_frac(0.90)
  test  <- dplyr::anti_join(d2, train, by = 'id')
  id_save[[i]] <- test$id
}




n <- nrow(d)
d$trt <- d$trt -1

unique(d$celltype)
sq <- as.numeric(d$celltype == 'squamous')
small <- as.numeric(d$celltype == 'smallcell')
large <- as.numeric(d$celltype == 'large')




X <- cbind(d$karno, d$prior, d$age, d$diagtime, d$trt, sq, small, large)
varname <- c('karno', 'prior','age', 'diagtime', 'treatment', 'sq', 'small', 'large')

colnames(X) <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8')
X[, 1] <- X[, 1] / 10
X[, 2] <- X[, 2] / 100
X[, 3] <- X[, 3] / 100
X[, 4] <- X[, 4] / 100
p <- ncol(X)

par(mar = c(4, 4, .8, 1))

dataframeX <- data.frame(X)
TMdataT <- cbind(d$time, d$status, X)
TMdataT <- data.frame(TMdataT)
colnames(TMdataT) <- c('time', 'delta', colnames(dataframeX))
TMdataT$time <- TMdataT$time/30


##### spBayesSurv setting 
nburn=4000; nsave=6000; nskip=15; niter = nburn+nsave;
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(maxL=15, a0=1, b0=1, M=10, q=.9);

### prediction evaluation
CindexBuLTM <- numeric(10)
CindexPO <- numeric(10)
CindexTM <- numeric(10)


IBSBuLTM <- numeric(10)
IBSPO <- numeric(10)
IBSTM <- numeric(10)

MAEBuLTM <- numeric(10)
MAEPO <- numeric(10)
MAETM <- numeric(10)

grids <- seq(.1, max(TMdataT$time), by = .1)

len <- length(grids)

# TMdataT$id <- 1:nrow(TMdataT)


for(j in 1:10){

# train <- TMdataT %>% dplyr::sample_frac(0.90)
# test  <- dplyr::anti_join(TMdataT, train, by = 'id')
test <- TMdataT[id_save[[j]], ]
train <- TMdataT[-id_save[[j]], ]

fit <- BuLTMfit(train$time, train[, -c(1, 2)], delta = train$delta, K = 15, iter = 1500, alpha = 1, eta = .1, probseries = seq(0, 1, by = .2))

S_train <- Surv(train$time, train$delta, type = "right")

fitPO <- survregbayes(S_train ~ x1+x2+x3+x4+x5+x6+x7+x8, data = train, survmodel = "PO", prior = prior, mcmc = mcmc,
                      dist = "loglogistic")

fitTM <- TransModel(S_train ~ x1+x2+x3+x4+x5+x6+x7+x8, data = train, r = 1)



BuLTM_med <- numeric(nrow(test))
PO_med <- numeric(nrow(test))
TM_med <- numeric(nrow(test))

sp_mat_BuLTM <- matrix(0, nrow(test), len)
sp_mat_PO <- matrix(0, nrow(test), len)
sp_mat_TM <- matrix(0, nrow(test), len)


for(i in 1:nrow(test)){
  predBuLTM <- predict(fit, X = as.numeric(test[i, 3:10]), grids = grids)
  BuLTM_med[i] <- grids[which.min(abs(predBuLTM$estimate - 0.5))]
  sp_mat_BuLTM[i, ] <- predBuLTM$estimate
  predPO <- plot.survregbayes(fitPO, tgrid = grids, xnewdata = test[i, 3:10], PLOT = F)
  PO_med[i] <- grids[which.min(abs(predPO$Shat - 0.5))]
  sp_mat_PO[i, ] <- predPO$Shat
  predTM <- predict(fitTM, new.time = grids, newdata = test[i, 3:10])
  TM_med[i] <- grids[which.min(abs(predTM$survival - 0.5))]
  sp_mat_TM[i, ] <- predTM$survival
  
  
  
  cat(i, '\n')
}


CindexBuLTM[j] <- Cindex(Surv(test$time, test$delta), BuLTM_med)
CindexPO[j] <- Cindex(Surv(test$time, test$delta), PO_med)
CindexTM[j] <- Cindex(Surv(test$time, test$delta), TM_med)


 IBSBuLTM[j] <- IBS(Surv(test$time, test$delta), sp_mat_BuLTM, grids)
 IBSPO[j] <- IBS(Surv(test$time, test$delta), sp_mat_PO, grids)
 IBSTM[j] <- IBS(Surv(test$time, test$delta), sp_mat_TM, grids)



MAEBuLTM[j] <- mean(abs(BuLTM_med - test$time)* test$delta)
MAEPO[j] <- mean(abs(PO_med - test$time)* test$delta)
MAETM[j] <- mean(abs(TM_med - test$time)* test$delta)

}


data_BS = data.frame('IBS' = c(IBSBuLTM, IBSPO, IBSTM),
                     'model' = c(rep('BuLTM', 10), rep('spBayesSurv', 10), rep('TransModel', 10))) 

ggplot(data_BS, aes(x = model, y = IBS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF", '#99FFCC'))+
  theme(legend.position="none") +
  theme(axis.title.x = element_blank())


data_CI = data.frame('Cindex' = c(CindexBuLTM, CindexPO, CindexTM),
                     'model' = c(rep('BuLTM', 10), rep('spBayesSurv', 10), rep('TransModel', 10))) 

ggplot(data_CI, aes(x = model, y = Cindex, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF", '#99FFCC')) +
  theme(legend.position="none") +
  theme(axis.title.x = element_blank())


data_MAE = data.frame('MAE' = c(MAEBuLTM, MAEPO, MAETM),
                      'model' = c(rep('BuLTM', 10), rep('spBayesSurv', 10), rep('TransModel', 10))) 

ggplot(data_MAE, aes(x = model, y = MAE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF", '#99FFCC')) +
  theme(legend.position="none") +
  theme(axis.title.x = element_blank())

# data_CI <- read.table('CindexAll.txt', header = T)
#data_MAE <-  read.table('MAEAll.txt', header = T)

t.test(CindexBuLTM, CindexPO, paired = T)
t.test(CindexBuLTM, CindexTM, paired = T)



#data_MAE %>% group_by(model) %>% summarise_at(vars(MAE), list(name = mean))

