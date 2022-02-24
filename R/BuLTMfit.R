# Save this file as `R/BuLTM.R`

#' Unified nonparametric Bayesian inference model with Stan
#'@import splines2
#'@import MASS
#'@import rstan
#'@import loo
#' @param X Numeric Matrix of input values
#' @param Yi Numeric Vector of output values
#' @param probseries Sequence from 0 to 1 or the step value of the sequence from 0 to 1
#' @param K Integer value of truncated DP
#' @param alpha Numeric value for mass parameter of DP
#' @param prior The prior for coefficient, defaut is 'Exp', exponential; otherwise truncated Gaussian
#' @param chains Argument passed to `rstan::sampling`.
#' @param iter Argument passed to `rstan::sampling` .
#' @param warmup Argument passed to `rstan::sampling` .
#' @param thin Argument passed to `rstan::sampling` .
#' @return An object of class `BuLTMfit` returned by `rstan::sampling`
#' @export
#
# library("bayesplot")
# library("ggplot2")
# library(rstan)
# library(splines2)
# library(survival)
# library(xtable)
# library(ggplot2)
# library(MASS)
# library(spBayesSurv)
# library(boot)
# library(HDInterval)
# library(distr)



######## model
BuLTMfit <- function(Yi,
                          X,
                          delta,
                          probseries,
                          K,
                          alpha,
                          prior = 'Exp',
                          eta = 1,
                          psi = 1,
                          v = 1,
                          chains = 4,
                          iter = 2500,
                          warmup = 500,
                          thin = 1
                    )
{
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  ###############  check input params
  if(!(delta==0 || delta==1)){
    stop("parameter 'delta' expects a boolean value")
  }

  probseries_len = length(probseries)
  if(probseries_len==1 && probseries>0 && probseries<1){
    probseries = seq(0, 1, probseries)
  }else if(probseries_len>1 && probseries[1]==0 && probseries[probseries_len]==1){
    probseries = probseries
  }else{
    stop("parameter 'probseries' expects a step value between 0 and 1, or a sequence from 0 to 1")
  }

  ###############  pre-process input params
  TMdataT <- as.data.frame(cbind(Yi,X, delta))
  ### knots -> for splines
  knots <- quantile(sort(Yi[which(delta == 1)]),
                    probs = probseries)
  knotdiff <- abs(ecdf(Yi)(knots) - probseries)
  if(max(knotdiff) > 0.05){
    knots <- sort(c(knots,quantile(Yi, probseries[which(knotdiff > 0.05)])))
  }
  ### t
  maxY <- max(Yi)
  t <- c(sort(Yi), seq(maxY, (maxY+0.05), .01))
  ### NN,location
  NN <- length(t)
  location <- numeric(length(Yi)) #The location of an observation in sorted time vector t
  for(i in 1:length(Yi)){
    location[i] <- which(t == Yi[i])[1]
  }
  ### splines
  Yspline <- t(iSpline(t, intercept = F, knots = knots, Boundary.knots = c(0, max(t)))) #spline of H(t)
  yspline <- t(iSpline(t, intercept = F, knots = knots, derivs = 1, Boundary.knots = c(0, max(t)))) #spline of h(t)
  num_basis <- nrow(Yspline)
  ### n,p
  n = nrow(X)
  p = ncol(X)
  ### beta_init0,w_init0
  beta_init0 <- rnorm(p)
  beta_init0 <- beta_init0 / sqrt(norm(beta_init0, '2'))
  w_init0 <- rbeta(K, 1, 1)

  ########### argument-data for rstan::sampling
  Model_data <- list(
    p=p,
    num_basis = num_basis,
    K = K,
    alpha = alpha,
    n=n,
    NN = NN,
    Yspline = Yspline,
    yspline = yspline,
    Y=Yi,
    location = location,
    delta=delta,
    X=X,
    eta = eta,
    psi = psi,
    v = v,
    beta_init0=beta_init0,
    w_init0=w_init0)

  ########### other params for rstan::sampling
  inits <- function()
    list(
      beta = beta_init0,

      w = w_init0
    )
  params <- c('beta', 'beta_trans', 'nu', 'theta', 'w', 'a_raw', 'DP_weights', 'log_lik')

  if(prior == 'Exp'){
    A <- rstan::sampling(stanmodels$trans_unknow_Exp,
                         data = Model_data,
                         chains = chains,
                         iter = iter,
                         warmup = warmup,
                         thin = thin)
  }
  else{
    A <- rstan::sampling(stanmodels$trans_unknow_Gaussian,
                         data = Model_data,
                         chains = chains,
                         iter = iter,
                         warmup = warmup,
                         thin = thin)
  }


  ########### process output
  chain <- extract(A)
  timing <- sum(get_elapsed_time(A))/nrow(get_elapsed_time(A))

  hat_beta_init <- apply(chain$beta, 2, mean)
  hat_beta <- apply(chain$beta_trans, 2, mean)

  MCMC_DIAG <- monitor(A, digits_summary = 4, print = F, warmup = warmup)

  log_lik <- extract_log_lik(A, merge_chains = F)
  r_eff <- relative_eff(exp(log_lik), cores = 4)
  Loo <- loo(log_lik, r_eff = r_eff, cores = 4)$estimate
  LPML <- Loo[1, 1]
  WAIC <- Loo[3, 1]

  z <- list()
  z$stan <- A
  z$chain <- chain
  z$hat_beta <- hat_beta
  z$MCMC_DIAG <- MCMC_DIAG
  z$t <- t
  z$K <- K
  z$knots <- knots
  z$LPML <- LPML
  z$WAIC <- WAIC
  z$time <- timing
  class(z) <- "BuLTMfit"
  return(z)
}



