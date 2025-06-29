#' Functions for sufficient informativeness criterion
#' @import rstan
#' @import ggplot2
#' @param object The BuLTM object
#' @param plot The logical variable for plot
#' @export

SuffitInf <- function(object, plot = TRUE){
  if(!inherits(object, "BuLTM"))
    stop("It is not a BuLTM object!")

  probseries <- object$probseries
  knots <- object$knots
  t <- object$t
  eta <- object$eta
  zeta <- object$zeta
  v <- object$v
  K <- object$K

  q0 <- min(probseries[which(probseries > 1- exp(-1))])
  t0 <- max(knots[names(knots) == paste0(as.character(q0*100), '%')])
  Yspline <- object$Yspline

  loc <- which.min(abs(t - t0))

  length <- object$length
  A1_chain <- object$chain

  c11 <- A1_chain$alpha[1:length, ] %*% Yspline[, loc]
  c12 <- A1_chain$alpha[(length+1):(2*length), ]  %*% Yspline[, loc]
  c13 <- A1_chain$alpha[(2*length+1):(3*length), ]  %*% Yspline[, loc]
  c14 <- A1_chain$alpha[(3*length+1):(4*length), ]  %*% Yspline[, loc]

  within_var <- mean(c(var(c11), var(c12), var(c13), var(c14)))


  knot.loc <- which(knots == t0)
  ess <- as.numeric(object$MCMC_DIAG['lp_','n_eff'])
  Rhat <- as.numeric(object$MCMC_DIAG['lp_','Rhat'])

  priorinfor_zeta <- function(zeta, knot.loc = 5, eta = .01, K=12, degree.spline = 4){
    1/(K*zeta + eta/knot.loc+degree.spline)^2 +
      1/(K*zeta + eta/(knot.loc+degree.spline))
  }


  # InfLevel <-  1/(K*zeta + eta/(knot.loc+4))^2 +
  #   1/(K*zeta + eta/(knot.loc+4))

  x_loc <- zeta
  y_loc <- priorinfor_zeta(x_loc, knot.loc = knot.loc,
                           eta = eta, K = K,
                           degree.spline = 4)
  if(!is.logical(plot))
    stop("plot must be TURE or FALSE!")

  candidate_zeta <- seq(.01, 1, by = .01)

  plot_data <- data.frame(
    candidate_zeta = candidate_zeta,
    priorinfor = sapply(candidate_zeta,
                        FUN = function(s){
                          priorinfor_zeta(
                            zeta = s, knot.loc = knot.loc,
                            eta = eta, K = K,
                            degree.spline = 4)})
  )

  p <- ggplot(plot_data, aes(x = candidate_zeta, y = priorinfor)) +
    geom_line(linewidth = 1.5) +
    geom_hline(yintercept = within_var, color = "red") +
    geom_vline(xintercept = x_loc, color = "blue", linetype = "dashed") +
    annotate("point",x = x_loc, y = y_loc, color = "green", size = 3) +
    annotate("text",
             x = x_loc,
             y = y_loc,
             label = sprintf("(%.2f, %.3f)", x_loc, y_loc),
             hjust = -0.5, vjust = -0.5, color = "black") +
    xlab(expression(zeta)) +
    ylab("Inverse of prior information level")
  if (plot) {
    print(p)
  }

  return(list(
    within_var = within_var,
    Inv_InfLevel = y_loc,
    knot.loc = knot.loc,
    ess = ess,
    Rhat = Rhat,
    y_loc = y_loc
  ))
}
