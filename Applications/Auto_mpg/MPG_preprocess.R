rm(list = ls())
set.seed(1234)

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


genFourier <- function(x, K = 10){
  Phi <- matrix(0, length(x), K)
  for(k in 1:K){
    if (k %% 2 == 1){
      l <- (k + 1) / 2
      Phi[, k] <-  sqrt(2) * sin((l) * pi * (2*x))
      
    }
    else{
      l <- k / 2
      Phi[, k] <-  sqrt(2) * cos(l * pi * (2*x))
    }
  }
  return(Phi)
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

X_Cont <- MPG[, 3:6]


non_lin_X <- c()

for(j in 1:ncol(X_Cont)){
  non_lin_X <- cbind(non_lin_X, genFourier(X_Cont[, j], K=8))
}



NonlinAug_Data <- cbind(MPG[, -(3:6)], non_lin_X)




name <- c()

name <- append(name, 'Y')

for(i in 1:(ncol(NonlinAug_Data) - 1)){
  name <- append(name, paste0('x', as.character(i)))
}


names(NonlinAug_Data) <- name


save(NonlinAug_Data, file = 'MPG_aug.RData')

