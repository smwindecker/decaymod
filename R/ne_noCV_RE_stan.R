#' Bayesian negative exponential function with Stan,
#' without cross validation, with random effects
#'
#' @export
#' @param mT Mass remaining of training data
#' @param m0 Initial mass of training data
#' @param time Vector of time values for training data
#' @param sp Numeric group identifier
#' @param J Number of groups
#' @param X Model matrix for training data#'
#' @param m0_sim Simulated log initial mass data
#' @param time_sim Simulated time data
#' @param sp_sim Numeric group identifier for sim dataset
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
ne_noCV_RE_stan <- function(mT, m0, time, sp, J, X, m0_sim, time_sim, sp_sim,...) {
  stan_data <- list(mT = mT,
                    m0 = m0,
                    time = time,
                    N = length(mT),
                    sp = sp,
                    J = J,
                    X = X,
                    P = ncol(X),
                    m0_sim = m0_sim,
                    time_sim = time_sim,
                    N_sim = length(time_sim),
                    sp_sim = sp_sim)
  out <- rstan::sampling(stanmodels$ne_noCV_RE, data = stan_data,
                         control = list(adapt_delta = 0.99, max_treedepth = 15), ...)
  return(out)
}
