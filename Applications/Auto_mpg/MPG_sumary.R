rm(list = ls())

library(ggplot2)
library(dplyr)


MPG_BuLTM_lin <- read.table('MPG_BuLTM_Lin.txt', header = T)
MPG_BuLTM_nonlin <- read.table('MPG_BuLTM.txt', header = T)
MPG_SeBR_lin <- read.table('MPG_SeBR_Lin.txt', header = T)
MPG_SeBR_nonlin <- read.table('MPG_SeBR.txt', header = T)
MPG_BCTM_lin <- read.table('MPG_BCTM_Lin.txt', header = T) 
MPG_mlt <- read.table('MPG_mlt_Lin.txt', header = T)
MPG_PTM <- read.csv('MPG_PTM.csv', header = T)


# boxplot(MPG_BuLTM_lin$CP, MPG_BuLTM_nonlin$CP, MPG_BCTM_lin$CP, MPG_SeBR_lin$CP, MPG_SeBR_nonlin$CP, MPG_mlt$CP, MPG_PTM$CP)
# 
# 
# abline(h = 0.95, col = 'red')
# 
# boxplot(MPG_BuLTM_lin$MAE, MPG_BuLTM_nonlin$MAE, MPG_BCTM_lin$MAE, MPG_SeBR_lin$MAE, MPG_SeBR_nonlin$MAE, MPG_mlt$MAE, MPG_PTM$MAE)
# 
# 


MAE <- c(MPG_BuLTM_lin$MAE, MPG_BuLTM_nonlin$MAE, MPG_SeBR_lin$MAE, MPG_SeBR_nonlin$MAE, 
         MPG_BCTM_lin$MAE, MPG_mlt$MAE, MPG_PTM$MAE)


CP <- c(MPG_BuLTM_lin$CP, MPG_BuLTM_nonlin$CP, MPG_SeBR_lin$CP, MPG_SeBR_nonlin$CP, 
        MPG_BCTM_lin$CP, MPG_mlt$CP, MPG_PTM$CP)


Model <- c(rep('BuLTM.lin', 10), rep('BuLTM.nonlin', 10), rep('SeBR.lin', 10), 
            rep('SeBR.nonlin', 10), rep('BCTM', 10), rep('mlt', 10), 
            rep('PTM', 10))


Dat <- data.frame(MAE = MAE, Model = factor(Model, levels = unique(Model)))

pic <- ggplot(Dat, aes(x=Model, y=MAE, fill = Model, alpha = .7)) +  geom_boxplot() + xlab("Model") +
  theme(legend.position="none") + xlab("") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "blue", linetype = "dashed")
pic



Dat <- data.frame(CP = CP, Model = factor(Model, levels = unique(Model)))

pic <- ggplot(Dat, aes(x=Model, y=CP, fill = Model, alpha = .7)) +  geom_boxplot() + xlab("Model") +
  theme(legend.position="none") + xlab("") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "blue", linetype = "dashed") +
 geom_hline(aes(yintercept=.95), colour="red", size = 1)

pic





t.test(MPG_BuLTM_lin$MAE, MPG_SeBR_lin$MAE, paired = T, alternative = 'less')
t.test(MPG_BuLTM_lin$MAE, MPG_mlt$MAE, paired = T)

