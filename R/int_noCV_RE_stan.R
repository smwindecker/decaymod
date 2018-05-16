#' Bayesian negative exponential plus intercept function with Stan,
#' without cross validation, with random effects
#'
#' @export
#' @param mT Mass remaining of training data
#' @param m0 Initial mass of training data
#' @param time Vector of time values for training data
#' @param sp Numeric group identifier
#' @param J Number of groups
#' @param X_k Model matrix for k for training data
#' @param X_i Model matrix for intercept for training data
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
int_noCV_RE_stan <- function(mT, m0, time, sp, J, X_k, X_i, ...) {
  stan_data <- list(mT = mT,
                    m0 = m0,
                    time = time,
                    N = length(mT),
                    sp = sp,
                    J = J,
                    X_k = X_k,
                    X_i = X_i,
                    P_k = ncol(X_k),
                    P_i = ncol(X_i))
  out <- rstan::sampling(stanmodels$int_noCV_RE, data = stan_data,
                         control = list(adapt_delta = 0.99, max_treedepth = 15), ...)
  return(out)
}
