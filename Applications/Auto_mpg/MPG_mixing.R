rm(list = ls())
library(dplyr)
library(splines2)
library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())

set.seed(1234)


load('MPG.RData')
Csigmoid <- function(x, C=5){
  if(C <= 0){stop('C is not positive')}
  else{
    C/(1+ exp(-x))
  }
}


priorinfor_zeta <- function(zeta, knot.loc = 5, eta = .01, K=12, degree.spline = 4){
  1/(K*zeta + eta/knot.loc+degree.spline)^2 + 
    1/(K*zeta + eta/(knot.loc+degree.spline))
}



d0 <- MPG


d0$id <- 1:nrow(d0)


  Y <- MPG$MPG

  X <- MPG[, -1]
  
  n <- nrow(MPG)
  
  delta <- rep(1, n)
 
  TY <- sapply(Y, FUN = function(x) Csigmoid(x))

  t <- c(sort(TY), seq(max(TY), (max(TY)+0.1), .01))
  NN <- length(t)
  
  location <- numeric(length(TY)) #The location of an observation in sorted time vector t 
  for(i in 1:length(TY)){
    location[i] <- which(t == TY[i])[1]
  }
  
  
  probseries <- seq(0, 1, .2)
  
  knots <- quantile(sort(TY[which(delta == 1)]), 
                    probs = probseries)
  
  
  Yspline <- t(iSpline(t, intercept = F, knots = knots, Boundary.knots = c(0, max(t)))) #spline ofn H(t) 
  yspline <- t(iSpline(t, intercept = F, knots = knots, derivs = 1, Boundary.knots = c(0, max(t)))) #spline of h(t)

  num_basis <- nrow(Yspline) #the dim of to be estimated coefficients
  
  
  
  
  a <- 1
  K <- 12
  p <- ncol(X)
  beta_init0 <- rnorm(p) #initial of \beta
  beta_init0 <- beta_init0 / norm(beta_init0, '2')
  w_init0 <- rbeta(K, 1, 1)
  
  params <- c('beta', 'psi', 'nu', 'w', 'alpha', 'beta_trans', 'DP_weights', 'log_lik')
  
  
  ###### determine the location
  t0 <- max(knots[names(knots) == '80%'])
  
  loc <- which.min(abs(t - t0))
  
    
  
  eta = .01
  zeta = .25
  v = 1
  
  
  Model_data <- list(
    p=p,
    num_basis = num_basis,
    K = K,
    a = a,
    n=n,
    NN = NN, 
    Yspline = Yspline,
    yspline = yspline,
    Y=TY,
    location = location,
    delta=delta,
    eta = eta, 
    zeta = zeta, 
    v = v, 
    X=X,
    beta_init0=beta_init0,
    w_init0=w_init0
  )
  inits <- function()
    list(
      beta = beta_init0,
      w = w_init0
    )
  
  A2 = stan(
    file = "trans_unknow_Exp.stan",
    data = Model_data,
    init = inits,
    #  algorithm =  "Fixed_param",
    pars = params,
    chains = 4,
    iter = 1000,
    warmup = 500,
    thin = 1, 
    seed = 1234
  )
  
  
  M2 <- monitor(A2, print = F)
  A2_chain <- extract(A2)

  
  c21 <- A2_chain$alpha[1:500, ] %*% Yspline[, loc]
  c22 <- A2_chain$alpha[501:1000, ] %*% Yspline[, loc]
  c23 <- A2_chain$alpha[1001:1500, ] %*% Yspline[, loc]
  c24 <- A2_chain$alpha[1501:2000, ] %*% Yspline[, loc]
  
  within_var <- mean(c(var(c21), var(c22), var(c23), var(c24))) 
  
  mcmc_trace(A2, pars = 'lp__')
  
  x_loc <- zeta
  y_loc <- priorinfor_zeta(x_loc)
  
  candidate_zeta <- seq(.1, 0.5, by = .05)
  
  plot_data <- data.frame(
    candidate_zeta = candidate_zeta,
    priorinfor = sapply(candidate_zeta, priorinfor_zeta)
  )
  
  
  
  ggplot(plot_data, aes(x = candidate_zeta, y = priorinfor)) +
    geom_line(linewidth = 1.5) +  # 注意：ggplot2 3.4.0+ 使用 linewidth 替代 lwd
    geom_hline(yintercept = within_var, color = "red") +
    geom_vline(xintercept = x_loc, color = "blue", linetype = "dashed") +
    geom_point(aes(x = x_loc, y = y_loc), color = "green", size = 3) +
    annotate("text", 
             x = x_loc, 
             y = y_loc,
             label = sprintf("(%.2f, %.3f)", x_loc, y_loc),
             hjust = -0.5, vjust = -0.5, color = "black") +
    xlab(expression(zeta)) +  
    ylab("Inverse of prior information level") 
  
  
  
  
  