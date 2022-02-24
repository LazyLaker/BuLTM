# Save this file as `R/BuLTM.R`

#' Unified nonparametric Bayesian inference model with Stan
#'@import splines2
#'@import MASS
#'@import rstan
#'@import loo
#' @param X Numeric Matrix of input values
#' @param Yi Numeric Vector of output values
#' @param num_knots number of knots
#' @param K Integer value of truncated DP
#' @param alpha Numeric value for mass parameter of DP
#' @param chains Argument passed to `rstan::sampling`.
#' @param iter Argument passed to `rstan::sampling` .
#' @param warmup Argument passed to `rstan::sampling` .
#' @param thin Argument passed to `rstan::sampling` .
#' @return An object of class `BuLTMfit2` returned by `rstan::sampling`
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
BuLTMfit2 <- function(Yi,
                     X,
                     delta,
                     num_knots=15,
                     K,
                     alpha,
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

  ###############  pre-process input params
  TMdataT <- as.data.frame(cbind(Yi,X, delta))
  ### t
  maxY <- max(Yi)
  t <- c(sort(Yi), seq(maxY, (maxY+0.05), .01))
  ### NN,location
  NN <- length(t)
  location <- numeric(length(Yi)) #The location of an observation in sorted time vector t
  for(i in 1:length(Yi)){
    location[i] <- which(t == Yi[i])[1]
  }
  ######knots for ispline
  knots <- seq(min(Yi), max(Yi), length.out = num_knots)

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
    beta_init0=beta_init0,
    w_init0=w_init0)

  ########### other params for rstan::sampling
  inits <- function()
    list(
      beta = beta_init0,

      w = w_init0
    )
  params <- c('beta', 'beta_trans', 'nu', 'theta', 'w', 'a_raw', 'DP_weights', 'log_lik')


  A <- rstan::sampling(stanmodels$trans_unknow_shrinkage,
                       data = Model_data,
                       chains = chains,
                       iter = iter,
                       warmup = warmup,
                       thin = thin)



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
  class(z) <- "BuLTMfit2"
  return(z)
}



