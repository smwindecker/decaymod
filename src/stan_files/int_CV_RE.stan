functions {
#include /functions/intercept_functions.stan
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
  int<lower=1> P_k;                   // num fixefs
  int<lower=1> P_i;                    // num fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P_k] X_k;         // design matrix for k effects
  matrix[1, P_k] X_k_test;         // design matrix for k effects
  matrix[J, P_i] X_i;           // design matrix for k effects
  matrix[1, P_i] X_i_test;
  int<lower=1, upper=J> sp[N];             // species id
}

parameters {
  vector[P_k] b_k;
  vector[P_i] b_i;
  vector[J] a_sp_k;                         // species intercepts on k
  vector[J] a_sp_i;                         // species intercepts on i
  real<lower=0> sigma_obs;                // observation sd
  real<lower=0> sigma_sp_k;                 // species sd
  real<lower=0> sigma_sp_i;
}

model {
  vector[J] k;
  vector[J] i;

  // priors
  a_sp_k ~ normal(0, sigma_sp_k);       // species random effects on k
  a_sp_i ~ normal(0, sigma_sp_i);         // species random effects on i
  sigma_sp_k ~ normal(0, 2);
  sigma_sp_i ~ normal(0, 2);
  sigma_obs ~ normal(0, 2);
  b_k ~ normal(0, 2);
  b_i ~ normal(0, 2);

  // likelihood
  k = param_re(b_k, X_k, P_k, J, a_sp_k);
  i = param_re(b_i, X_i, P_i, J, a_sp_i);
  mT ~ normal(m0 - k[sp] .* time + i[sp], sigma_obs);
}

generated quantities {

  real neg_loglik;
  vector[N_test] mT_pred;

  {
  real k_pred;
  real i_pred;

  k_pred = param_pred(b_k, X_k_test, P_k);
  i_pred = param_pred(b_i, X_i_test, P_i);
  mT_pred = intercept_pred_rng(N_test, m0_test, time_test, i_pred, k_pred, sigma_obs);
  neg_loglik = intercept_negloglik(mT_pred, m0_test, time_test, i_pred, k_pred, sigma_obs, N_test);
  }
}

