functions {
#include /functions/neg_exp_functions.stan
#include /functions/param_functions.stan
}

data {
  int<lower=1> N;                         // number of data points
  vector[N] mT;                             // logged mass at time T
  vector[N] m0;                             // logged initial mass
  vector<lower=0>[N] time;                  // time
  int<lower=1> P;                   // no fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P] X;         // design matrix for alpha effects
  int<lower=1, upper=J> sp[N];            // species id
}

parameters {
  vector[P] b;
  real<lower=0> sigma_obs;                // observation sd
}

model {
  vector[J] k;

  // priors
  sigma_obs ~ normal(0, 2);
  b ~ normal(0, 2);

  // likelihood
  k = param(b, X, P, J);
  mT ~ normal(m0 - k[sp] .* time, sigma_obs); // element-wise multiplication of vectors
}

generated quantities {

  vector[N] mT_fit;
  vector[J] k_fit;

  k_fit = param(b, X, P, J);
  mT_fit = negexp_fit_rng(N, m0, time, k_fit, sp, sigma_obs);

}
