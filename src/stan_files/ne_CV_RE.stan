functions {
#include /functions/neg_exp_functions.stan
#include /functions/param_functions.stan
}

data {
  int<lower=1> N;                         // number of data points
  int<lower=1> N_test;
  vector[N] mT;                             // logged mass at time T
  vector[N_test] mT_test;
  vector[N] m0;                             // logged initial mass
  vector[N_test] m0_test;
  vector<lower=0>[N] time;                  // time
  vector<lower=0>[N_test] time_test;
  int<lower=1> P;                   // no fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P] X;         // design matrix for alpha effects
  matrix[1, P] X_test;
  int<lower=1, upper=J> sp[N];            // species id
}

parameters {
  vector[P] b;
  vector[J] a;                         // species intercepts on alpha
  real<lower=0> sigma_obs;                // observation sd
  real<lower=0> sigma_sp;                 // species sp
}

model {
  vector[J] k;

  // priors
  a ~ normal(0, sigma_sp);       // species random effects on alpha
  sigma_sp ~ normal(0, 2);
  sigma_obs ~ normal(0, 2);
  b ~ normal(0, 2);

  // likelihood
  k = param_re(b, X, P, J, a);
  mT ~ normal(m0 - k[sp] .* time, sigma_obs);
}

generated quantities {

  real neg_loglik;
  vector[N_test] mT_pred;

  {
  real k_pred;

  k_pred = param_pred(b, X_test, P);
  mT_pred = negexp_pred_rng(N_test, m0_test, time_test, k_pred, sigma_obs);
  neg_loglik = negexp_negloglik(mT_test, m0_test, time_test, k_pred, sigma_obs, N_test);
  }
}
