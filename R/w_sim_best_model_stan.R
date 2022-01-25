#' Bayesian weibull function with Stan applying best model to simulated
#' data
#'
#' @export
#' @param mT Mass remaining of training data
#' @param m0 Initial mass of training data
#' @param time Vector of time values for training data
#' @param sp Numeric group identifier
#' @param J Number of groups
#' @param X_alpha Model matrix for alpha for training data
#' @param X_beta Model matrix for beta for training data
#' @param X_alpha_sim Simulated model matrix for alpha
#' @param X_beta_sim Simulated model matrix for beta
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
w_sim_best_model_stan <- function(mT, m0, time, sp, J, X_alpha,
                                  X_beta, X_alpha_sim, X_beta_sim, ...) {

  stan_data <- list(mT = mT,
                    m0 = m0,
                    time = time,
                    N = length(mT),
                    sp = sp,
                    J = J,
                    X_alpha = X_alpha,
                    X_beta = X_beta,
                    P_alpha = ncol(X_alpha),
                    P_beta = ncol(X_beta),
                    m0_sim = log(4100),
                    time_sim = 0.7,
                    J_sim = nrow(X_beta_sim),
                    X_alpha_sim = X_alpha_sim,
                    X_beta_sim = X_beta_sim,
                    sp_sim = 1:10000)
  out <- rstan::sampling(stanmodels$w_sim_best_model, data = stan_data,
                         control = list(adapt_delta = 0.99, max_treedepth = 15), ...)
  return(out)
}
