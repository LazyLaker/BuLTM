# Save this file as `R/BuLTM.R`

#' Bayesian unidentified Linear Transformation model
#'@import splines2
#'@import MASS
#'@import rstan
#'@import loo
#' @param X Numeric Matrix of input values.
#' @param Y Numeric Vector of output values.
#' @param probseries Sequence from 0 to 1 or the step value of the sequence from 0 to 1.
#' @param K Integer value of truncated DP.
#' @param a Numeric value for mass parameter of DP.
#' @param eta hyperparameter eta.
#' @param zeta hyperparameter zeta.
#' @param v hyperparameter v.
#' @param chains Argument passed to `rstan::sampling`.
#' @param iter Argument passed to `rstan::sampling` .
#' @param warmup Argument passed to `rstan::sampling` .
#' @param thin Argument passed to `rstan::sampling` .
#' @param seed The seed for sampling.
#' @param adapt_delta The adapt_delta in Stan.
#' @param max_treedepth The max_treedepth in Stan.
#' @return An object of class `BuLTM` returned by `rstan::sampling`
#' @export
#




######## model
BuLTM <- function(Y,
                     X,
                     delta,
                     probseries,
                     K,
                     a = 1,
                     eta = .01,
                     zeta = 1,
                     v = 1,
                     chains = 4,
                     iter = 1500,
                     warmup = 500,
                     thin = 1,
                     seed = 1234,
                     adapt_delta = 0.8,
                     max_treedepth = 9
)
{
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  ###############  check input params


  probseries_len = length(probseries)
  if(probseries_len==1 && probseries>0 && probseries<1){
    probseries = seq(0, 1, probseries)
  }else if(probseries_len>1 && probseries[1]==0 && probseries[probseries_len]==1){
    probseries = probseries
  }else{
    stop("parameter 'probseries' expects a step value between 0 and 1, or a sequence from 0 to 1")
  }

  ###############  pre-process input params
  TMdataT <- as.data.frame(cbind(Y,X, delta))
  ### knots -> for splines
  knots <- quantile(sort(Y[which(delta == 1)]),
                    probs = probseries)
  knotdiff <- abs(ecdf(Y)(knots) - probseries)
  if(max(knotdiff) > 0.05){
    knots <- sort(c(knots,quantile(Y, probseries[which(knotdiff > 0.05)])))
  }
  ### t
  maxY <- max(Y)
  t <- c(sort(Y), seq(maxY, (maxY+0.05), .01))
  ### NN,location
  NN <- length(t)
  location <- numeric(length(Y)) #The location of an observation in sorted time vector t
  for(i in 1:length(Y)){
    location[i] <- which(t == Y[i])[1]
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
    a = a,
    n=n,
    NN = NN,
    Yspline = Yspline,
    yspline = yspline,
    Y=Y,
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
  ########### other params for rstan::sampling
  params <- c('beta', 'beta_trans', 'nu', 'theta', 'w', 'a_raw', 'DP_weights', 'log_lik')


  A <- rstan::sampling(stanmodels$trans_unknow_Exp,
                         data = Model_data,
                         chains = chains,
                         iter = iter,
                         warmup = warmup,
                         thin = thin, seed = seed,
                         control = list(adapt_delta = adapt_delta,
                                        max_treedepth = max_treedepth)
                         )


  ########### process output
  chain <- extract(A)
  timing <- sum(get_elapsed_time(A))/nrow(get_elapsed_time(A))

  hat_beta_init <- apply(chain$beta, 2, mean)
  hat_beta <- apply(chain$beta_trans, 2, mean)

  MCMC_DIAG <- monitor(A, digits_summary = 4, print = F, warmup = warmup)

  log_lik <- extract_log_lik(A, merge_chains = F)
  r_eff <- relative_eff(exp(log_lik), cores = 4)
  Loo <- loo(log_lik, r_eff = r_eff, cores = 4)$estimate
  WAIC <- Loo[3, 1]


  z <- list()
  z$stan <- A
  z$chain <- chain
  z$hat_beta <- hat_beta
  z$MCMC_DIAG <- MCMC_DIAG
  z$t <- t
  z$K <- K
  z$eta <- eta
  z$zeta <- zeta
  z$v <- v
  z$Yspline <- Yspline
  z$knots <- knots
  z$WAIC <- WAIC
  z$probseries <- probseries
  z$time <- timing
  z$length <- iter - warmup
  class(z) <- "BuLTM"
  return(z)
}



