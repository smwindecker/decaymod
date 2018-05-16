functions {
#include /functions/intercept_functions.stan
#include /functions/param_functions.stan
}

data {
  int<lower=1> N;                         // number of data points
  vector[N] mT;                             // logged mass at time T
  vector[N] m0;                             // logged initial mass
  vector<lower=0>[N] time;                  // time
  int<lower=1> P_k;                   // num fixefs
  int<lower=1> P_i;                    // num fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P_k] X_k;         // design matrix for k effects
  matrix[J, P_i] X_i;           // design matrix for k effects
  int<lower=1, upper=J> sp[N];            // species id
}

parameters {
  vector[P_k] b_k;
  vector[P_i] b_i;
  real<lower=0> sigma_obs;                // observation sd
}

model {
  vector[J] k;
  vector[J] i;

  // priors
  sigma_obs ~ normal(0, 2);
  b_k ~ normal(0, 2);
  b_i ~ normal(0, 2);

  // likelihood
  k = param(b_k, X_k, P_k, J);
  i = param(b_i, X_i, P_i, J);
  mT ~ normal(m0 - k[sp] .* time + i[sp], sigma_obs);
}

generated quantities {

  vector[N] mT_fit;

  {
  vector[J] k_fit;
  vector[J] i_fit;

  k_fit = param(b_k, X_k, P_k, J);
  i_fit = param(b_i, X_i, P_i, J);
  mT_fit = intercept_fit_rng(N, m0, time, i_fit, k_fit, sp, sigma_obs);
  }
}
