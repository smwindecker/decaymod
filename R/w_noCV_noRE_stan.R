#' Bayesian weibull function with Stan,
#' without cross validation, without random effects
#'
#' @export
#' @param mT Mass remaining of training data
#' @param m0 Initial mass of training data
#' @param time Vector of time values for training data
#' @param sp Numeric group identifier
#' @param J Number of groups
#' @param X_alpha Model matrix for alpha for training data
#' @param X_beta Model matrix for beta for training data
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
w_noCV_noRE_stan <- function(mT, m0, time, sp, J, X_alpha, X_beta, ...) {
  stan_data <- list(mT = mT,
                    m0 = m0,
                    time = time,
                    N = length(mT),
                    sp = sp,
                    J = J,
                    X_alpha = X_alpha,
                    X_beta = X_beta,
                    P_alpha = ncol(X_alpha),
                    P_beta = ncol(X_beta))
  out <- rstan::sampling(stanmodels$w_noCV_noRE, data = stan_data,
                         control = list(adapt_delta = 0.99, max_treedepth = 15), ...)
  return(out)
}
