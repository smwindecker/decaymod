functions {
#include /functions/weibull_functions.stan
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
  int<lower=1> P_alpha;                   // num fixefs
  int<lower=1> P_beta;                    // num fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P_alpha] X_alpha;         // design matrix for alpha effects
  matrix[1, P_alpha] X_alpha_test;         // design matrix for alpha effects
  matrix[J, P_beta] X_beta;           // design matrix for alpha effects
  matrix[1, P_beta] X_beta_test;
  int<lower=1, upper=J> sp[N];             // species id
}

parameters {
  vector[P_alpha] b_alpha;
  vector[P_beta] b_beta;
  vector[J] a_sp_alpha;                         // species intercepts on alpha
  vector[J] a_sp_beta;                         // species intercepts on beta
  real<lower=0> sigma_obs;                // observation sd
  real<lower=0> sigma_sp_alpha;                 // species sd
  real<lower=0> sigma_sp_beta;
}

model {
  vector[J] alpha;
  vector[J] beta;
  vector[N] mu;

  // priors
  a_sp_alpha ~ normal(0, sigma_sp_alpha);       // species random effects on alpha
  a_sp_beta ~ normal(0, sigma_sp_beta);         // species random effects on beta
  sigma_sp_alpha ~ normal(0, 2);
  sigma_sp_beta ~ normal(0, 2);
  sigma_obs ~ normal(0, 2);
  b_alpha ~ normal(0, 2);
  b_beta ~ normal(0, 2);

  // likelihood
  alpha = param_re(b_alpha, X_alpha, P_alpha, J, a_sp_alpha);
  beta = param_re(b_beta, X_beta, P_beta, J, a_sp_beta);
  mu = weibull(N, m0, time, beta, alpha, sp, sigma_obs);
  mT ~ normal(mu, sigma_obs);

}

generated quantities {

  real neg_loglik;
  vector[N_test] mT_pred;

  {
  real alpha_pred;
  real beta_pred;

  alpha_pred = param_pred(b_alpha, X_alpha_test, P_alpha);
  beta_pred = param_pred(b_beta, X_beta_test, P_beta);
  mT_pred = weibull_pred_rng(N_test, m0_test, time_test, beta_pred, alpha_pred, sigma_obs);
  neg_loglik = weibull_negloglik(mT_pred, m0_test, time_test, beta_pred, alpha_pred, sigma_obs, N_test);
  }
}

