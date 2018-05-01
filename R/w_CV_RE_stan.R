#' Bayesian weibull function with Stan,
#' including cross validation, with random effects
#'
#' @export
#' @param mT Mass remaining of training data
#' @param mT_test Mass remaining of test data
#' @param m0 Initial mass of training data
#' @param m0_test Initial mass of test data
#' @param time Vector of time values for training data
#' @param time_test Vector of time values for test data
#' @param sp Numeric group identifier
#' @param J Number of groups
#' @param X_alpha Model matrix for alpha for training data
#' @param X_alpha_test model matrix for alpha for test data
#' @param X_beta Model matrix for beta for training data
#' @param X_beta_test model matrix for beta for test data
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
w_CV_RE_stan <- function(mT, mT_test, m0, m0_test, time, time_test, sp, J,
                         X_alpha, X_alpha_test, X_beta, X_beta_test, ...) {
  stan_data <- list(mT = mT,
                    mT_test = mT_test,
                    m0 = m0,
                    m0_test = m0_test,
                    time = time,
                    time_test = time_test,
                    N = length(mT),
                    N_test = length(mT_test),
                    sp = sp,
                    J = J,
                    X_alpha = X_alpha,
                    X_alpha_test = X_alpha_test,
                    X_beta = X_beta,
                    X_beta_test = X_beta_test,
                    P_alpha = ncol(X_alpha),
                    P_beta = ncol(X_beta))
  out <- rstan::sampling(stanmodels$w_CV_RE, data = stan_data,
                         control = list(adapt_delta = 0.99, max_treedepth = 15), ...)
  return(out)
}
