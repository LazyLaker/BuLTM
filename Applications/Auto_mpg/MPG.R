rm(list = ls())
set.seed(1234)
library(dplyr)
library(caret)

MPG <- read.table('auto-mpg.data', dec = ',', header = T)

name <- c('MPG', 'cylinders', 'displacement', 'horsrpower', 
          'weight', 'acceleration', 'modelyear', 'origin', 'name')

names(MPG) <- name

unif_trans <- function(x){
  maxx <- max(x)
  minx <- min(x)
  res <- sapply(x, FUN = function(z) (z-minx)/(maxx - minx))
  return(res)
}


MPG <- MPG[, -ncol(MPG)]  


MPG$MPG <- as.numeric(MPG$MPG)  
MPG$displacement <- as.numeric(MPG$displacement)
MPG$horsrpower <- as.numeric(MPG$horsrpower)
MPG$weight <- as.numeric(MPG$weight)
MPG$acceleration <- as.numeric(MPG$acceleration)

MPG <- na.omit(MPG)

n <- nrow(MPG)


MPG$MPG <- as.numeric(scale(MPG$MPG))

MPG$displacement <- unif_trans(MPG$displacement)
MPG$horsrpower <- unif_trans(MPG$horsrpower)

MPG$weight <- unif_trans(MPG$weight)
MPG$acceleration <- unif_trans(MPG$acceleration)


save(MPG, file = 'MPG.RData')


write.csv(MPG, file = 'MPG.csv', row.names = F)



######## seperate test set 


id_save <- list()

d0 <- MPG


d0$id <- 1:nrow(d0)
m <- 10
for(i in 1:m){
  train <- d0 %>% dplyr::sample_frac(0.90)
  test  <- dplyr::anti_join(d0, train, by = 'id')
  id_save[[i]] <- test$id
}

test_index <- do.call(rbind, id_save)


write.csv(data.frame(test_index), file = 'test_index.csv', row.names = F)


