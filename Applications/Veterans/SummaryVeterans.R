
library(ggplot2)

library(dplyr)
library(tidyverse)
A <- read.table('VeteransCV.txt',  header = T)

n <- nrow(A)


Cindex <- c(A$CindexBuLTM, A$CindexPO, A$CindexTM, A$Cindmlt)

Model <- c(rep('BuLTM', n), rep('spBayesSurv', n), rep('TransModel', n), rep('mlt', n))

Dat <- data.frame(Cindex = Cindex, Model = factor(Model, levels = unique(Model)))

pic <- ggplot(Dat, aes(x=Model, y=Cindex, fill = Model, alpha = 0.7)) +  geom_boxplot() + xlab("Model") +
  theme(legend.position="none") + xlab("") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "red", linetype = "dashed")
pic



t.test(A$CindexBuLTM, A$Cindmlt, paired = T)


IBS <- c(A$IBSBuLTM, A$IBSPO, A$IBSTM, A$IBSmlt)

IBS[which(IBS > 1)] <- NA

IBS <- na.omit(IBS)

n <- length(IBS)/4

Model <- c(rep('BuLTM', n), rep('spBayesSurv', n), rep('TransModel', n), rep('mlt', n))

Dat <- data.frame(IBS = IBS, Model = factor(Model, levels = unique(Model)))

pic <- ggplot(Dat, aes(x=Model, y=IBS, fill = Model, alpha = 0.7)) +  geom_boxplot() + xlab("Model") +
  theme(legend.position="none") + xlab("") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "red", linetype = "dashed")
pic

IBS_BuLTM <- Dat %>% filter(Model == 'BuLTM')
IBS_spBayesSurv <- Dat %>% filter(Model == 'spBayesSurv')
IBS_tm <- Dat %>% filter(Model == 'TransModel')
IBS_mlt <- Dat %>% filter(Model == 'mlt')




t.test(IBS_BuLTM$IBS, IBS_spBayesSurv$IBS, paired = T)
t.test(IBS_BuLTM$IBS, IBS_tm$IBS, paired = T)

t.test(IBS_BuLTM$IBS, IBS_mlt$IBS, paired = T)
