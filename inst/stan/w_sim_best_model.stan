functions {
#include /functions/weibull_functions.stan
#include /functions/param_functions.stan
}

data {
  int<lower=1> N;                         // number of data points
  vector[N] mT;                             // logged mass at time T
  vector[N] m0;                             // logged initial mass
  vector<lower=0>[N] time;                  // time
  int<lower=1> P_alpha;                   // num fixefs
  int<lower=1> P_beta;                    // num fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P_alpha] X_alpha;         // design matrix for alpha effects
  matrix[J, P_beta] X_beta;           // design matrix for alpha effects
  int<lower=1, upper=J> sp[N];            // species id

  int<lower=1> N_sim;                         // number of data points
  real<lower=0> m0_sim;
  real<lower=0> time_sim;                   // time
  int<lower=1> P_alpha_sim;                   // num fixefs
  int<lower=1> P_beta_sim;                    // num fixefs
  int<lower=1> J_sim;                         // number of species
  matrix[J_sim, P_alpha_sim] X_alpha_sim;         // design matrix for alpha effects
  matrix[J_sim, P_beta_sim] X_beta_sim;           // design matrix for alpha effects
  int<lower=1, upper=J_sim> sp_sim[N_sim];            // species id

}

parameters {
  vector[P_alpha] b_alpha;
  vector[P_beta] b_beta;
  real<lower=0> sigma_obs;                // observation sd
}

model {
  vector[J] alpha;
  vector[J] beta;
  vector[N] mu;

  // priors
  sigma_obs ~ normal(0, 2);
  b_alpha ~ normal(0, 2);
  b_beta ~ normal(0, 2);

  // likelihood
  alpha = param(b_alpha, X_alpha, P_alpha, J);
  beta = param(b_beta, X_beta, P_beta, J);
  mu = weibull(N, m0, time, beta, alpha, sp, sigma_obs);
  mT ~ normal(mu, sigma_obs);
}

generated quantities {

  vector[N_sim] mT_sim;
  vector[J_sim] alpha_sim;
  vector[J_sim] beta_sim;

  alpha_sim = param(b_alpha, X_alpha_sim, P_alpha_sim, J_sim);
  beta_sim = param(b_beta, X_beta_sim, P_beta_sim, J_sim);
  mT_sim = weibull_best_sim_rng(N_sim, m0_sim, time_sim, beta_sim, alpha_sim, sp_sim);

}
