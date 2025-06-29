#' Predict conditional survival functions for Unified nonparametric Bayesian inference model
#'@import MASS
#'@import rstan
#' @param object The BuLTM object
#' @param X vector of covariates
#' @param grids The time grids to be predicted
#' @return A list including the estimated curve and 95% intervals
#' @export



predict.BuLTM <- function(object, X=NULL, grids = NULL)
{
  if(!inherits(object, "BuLTM"))
    warning("calling predict.BuLTM(<fake-BuLTM-object>) ...")
  if(is.numeric(X) == 0) stop('Only numeric covariates available')
  p <- ncol(object$chain$beta)
  if(length(X) != p) stop('Unavailable dimensions of X')
  if(is.null(X)) {X <- rep(0, p)}
  K <- object$K
  t <- object$t
  knots <- object$knots
  chain_beta <- object$chain$beta

  alpha <- object$chain$alpha
  DP_weights <- object$chain$DP_weights
  psi <- object$chain$psi
  nu <- object$chain$nu
  if(is.null(grids)){
    grids <- seq(.1, max(t), by = .01)
  }
  gridsspline <- t(iSpline(grids, intercept = F, knots = knots, Boundary.knots = c(0, max(grids)+0.01)))

  sample_path <- matrix(0, nrow(alpha), length(grids))

  temp_q <- matrix(rep(exp(-chain_beta %*% X), each = ncol(gridsspline)), ncol = ncol(gridsspline), byrow = T) * (alpha %*% gridsspline)
  for (i in 1:nrow(alpha)) {
    res <- rep(0, length(grids))
    for(k in 1:K){
      res <- res + DP_weights[i, k]*pweibull(matrix(temp_q[i,]), psi[i, k], nu[i, k])
    }
    sample_path[i, ] <- 1-res
  }

  pathmean <- apply(sample_path, 2, mean)
  pathmean_high <- apply(sample_path, 2, function(x) quantile(x, .975))
  pathmean_low <- apply(sample_path, 2, function(x) quantile(x, .025))
  res <- list(
    estimate = pathmean,
    high = pathmean_high,
    low = pathmean_low
  )
  return(res)
}
