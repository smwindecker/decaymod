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

  vector[N] mT_fit;
  vector[J] alpha_fit;
  vector[J] beta_fit;

  alpha_fit = param(b_alpha, X_alpha, P_alpha, J);
  beta_fit = param(b_beta, X_beta, P_beta, J);
  mT_fit = weibull_fit_rng(N, m0, time, beta_fit, alpha_fit, sp, sigma_obs);

}
