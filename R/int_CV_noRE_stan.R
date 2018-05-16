#' Bayesian negative exponential plus intercept function with Stan,
#' including cross validation, without random effects
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
#' @param X_k Model matrix for k for training data
#' @param X_k_test model matrix for k for test data
#' @param X_i Model matrix for intercept for training data
#' @param X_i_test model matrix for intercept for test data
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
int_CV_noRE_stan <- function(mT, mT_test, m0, m0_test, time, time_test, sp, J,
                           X_k, X_k_test, X_i, X_i_test, ...) {
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
                    X_k = X_k,
                    X_k_test = X_k_test,
                    X_i = X_i,
                    X_i_test = X_i_test,
                    P_k = ncol(X_k),
                    P_i = ncol(X_i))
  out <- rstan::sampling(stanmodels$int_CV_noRE, data = stan_data,
                         control = list(adapt_delta = 0.99, max_treedepth = 15), ...)
  return(out)
}
