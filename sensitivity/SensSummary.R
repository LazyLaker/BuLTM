library(ggplot2)

##### Box Gauss

n <- 100
Gauss_10 <- read.table('/Users/zhongchong/OneDrive - The Hong Kong Polytechnic University/Biometrika_BuLTM/sensitivity/sensitivity/result-0.1/May2025_Gauss.txt', header = T)
Gauss_4 <- read.table('/Users/zhongchong/OneDrive - The Hong Kong Polytechnic University/Biometrika_BuLTM/sensitivity/sensitivity/result-0.25/May2025_Gauss.txt', header = T)
Gauss_5 <- read.table('/Users/zhongchong/OneDrive - The Hong Kong Polytechnic University/Biometrika_BuLTM/simluation/CTM_2/Apr232025_Gauss.txt', header = T)


RIMSE <- c(Gauss_10$IMSEBuLTM, Gauss_5$IMSEBuLTM, Gauss_4$IMSEBuLTM)
MAE <- c(Gauss_10$MAEBuLTM, Gauss_5$MAEBuLTM, Gauss_4$MAEBuLTM)

t.test(Gauss_4$IMSEBuLTM, Gauss_5$IMSEBuLTM, paired = T)


KnotNumber <- c(rep(10, n), rep(5, n), rep(4, n))

Dat <- data.frame(RIMSE = RIMSE, Model = factor(KnotNumber, levels = unique(KnotNumber)))


pic <- ggplot(Dat, aes(x=factor(KnotNumber), y=RIMSE, fill = Model, alpha = .7)) +  geom_boxplot() + xlab("Knot Number") +
  theme(legend.position="none") + xlab("Number of knots") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "red", linetype = "dashed")
pic


Dat <- data.frame(MAE = MAE, Model = factor(KnotNumber, levels = unique(KnotNumber)))


pic <- ggplot(Dat, aes(x=factor(KnotNumber), y=MAE, fill = Model, alpha = .7)) +  geom_boxplot() + xlab("Knot Number") +
  theme(legend.position="none") + xlab("Number of knots") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "red", linetype = "dashed")
pic


##### box mix


n <- 100
Gauss_10 <- read.table('/Users/zhongchong/OneDrive - The Hong Kong Polytechnic University/Biometrika_BuLTM/sensitivity/sensitivity/result-0.1/May2025_GaussMix.txt', header = T)
Gauss_4 <- read.table('/Users/zhongchong/OneDrive - The Hong Kong Polytechnic University/Biometrika_BuLTM/sensitivity/sensitivity/result-0.25/May2025_GaussMix.txt', header = T)
Gauss_5 <- read.table('/Users/zhongchong/OneDrive - The Hong Kong Polytechnic University/Biometrika_BuLTM/simluation/CTM_2/Apr152025_GaussMix.txt', header = T)




RIMSE <- c(Gauss_10$IMSEBuLTM, Gauss_5$IMSEBuLTM, Gauss_4$IMSEBuLTM)
MAE <- c(Gauss_10$MAEBuLTM, Gauss_5$MAEBuLTM, Gauss_4$MAEBuLTM)


KnotNumber <- c(rep(10, n), rep(5, n), rep(4, n))

Dat <- data.frame(RIMSE = RIMSE, Model = factor(KnotNumber, levels = unique(KnotNumber)))


pic <- ggplot(Dat, aes(x=factor(KnotNumber), y=RIMSE, fill = Model, alpha = .7)) +  geom_boxplot() + xlab("Knot Number") +
  theme(legend.position="none") + xlab("Number of knots") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "red", linetype = "dashed")
pic


t.test(Gauss_10$IMSEBuLTM, Gauss_4$IMSEBuLTM, paired = T)


Dat <- data.frame(MAE = MAE, Model = factor(KnotNumber, levels = unique(KnotNumber)))


pic <- ggplot(Dat, aes(x=factor(KnotNumber), y=MAE, fill = Model, alpha = .7)) +  geom_boxplot() + xlab("Knot Number") +
  theme(legend.position="none") + xlab("Number of knots") + stat_summary(fun = "mean",
                                                          geom = "crossbar", size = .2,
                                                          color = "red", linetype = "dashed")
pic







