functions {
#include /functions/weibull_functions.stan
#include /functions/param_functions.stan
}

data {
  int<lower=1> N;                         // number of data points
  vector[N] mT;                             // logged mass at time T
  vector[N] m0;                             // logged initial mass
  vector<lower=0>[N] time;                   // time
  int<lower=1> P_alpha;                   // no fixefs
  int<lower=1> P_beta;                    // no fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P_alpha] X_alpha;             // design matrix for alpha effects
  matrix[J, P_beta] X_beta;               // design matrix for alpha effects
  int<lower=1, upper=J> sp[N];            // species id

  int<lower=1> N_sim;                         // number of data points
  vector[N_sim] m0_sim;
  vector<lower=0>[N_sim] time_sim;                   // time
  int<lower=1, upper=J> sp_sim[N_sim];            // species id
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

  vector[N_sim] mT_sim;
  vector[J] alpha_fit;
  vector[J] beta_fit;

  alpha_fit = param_re(b_alpha, X_alpha, P_alpha, J, a_sp_alpha);
  beta_fit = param_re(b_beta, X_beta, P_beta, J, a_sp_beta);
  mT_sim = weibull_sim_rng(N_sim, m0_sim, time_sim, beta_fit, alpha_fit, sp_sim);

}
