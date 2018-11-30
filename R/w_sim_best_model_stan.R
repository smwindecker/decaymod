#' Bayesian weibull function with Stan,
#' without cross validation, without random effects
#' with simulation for traits
#'
#' @export
#' @param mT Log mass remaining of training data
#' @param m0 Log initial mass of training data
#' @param time Vector of time values for training data
#' @param sp Numeric group identifier
#' @param J Number of groups
#' @param X_alpha Model matrix for alpha for training data
#' @param X_beta Model matrix for beta for training data
#' @param m0_sim Simulated log initial mass data
#' @param time_sim Simulated time data
#' @param J_sim Number of simulated groups
#' @param X_alpha_sim Model maxtrix for simulated alpha
#' @param X_beta_sim Model matrix for simulated beta
#' @param sp_sim Numeric group identifier for sim dataset
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
w_sim_best_model_stan <- function (mT, m0, time, sp, J, X_alpha, X_beta, m0_sim,
                                  time_sim, J_sim, X_alpha_sim, X_beta_sim, sp_sim, ...) {
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
                    m0_sim = m0_sim,
                    time_sim = time_sim,
                    sp_sim = sp_sim,
                    J_sim = J_sim,
                    N_sim = J_sim,
                    X_alpha_sim = X_alpha_sim,
                    X_beta_sim = X_beta_sim,
                    P_alpha_sim = ncol(X_alpha_sim),
                    P_beta_sim = ncol(X_beta_sim))

  out <- rstan::sampling(stanmodels$w_sim_best_model, data = stan_data,
                         control = list(adapt_delta = 0.99, max_treedepth = 15), ...)
  return(out)
}
