
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
  int<lower=1, upper=J> sp[N];            // species id
}

parameters {
  vector[P_alpha] b_alpha;
  vector[P_beta] b_beta;
  real<lower=0> sigma_obs;                // observation sd
}

model {

  real alpha[J];
  real alpha_ln[J];
  real beta[J];
  real beta_ln[J];
  vector[N] mu;

  // priors
  sigma_obs ~ normal(0, 2);
  b_alpha ~ normal(0, 2);
  b_beta ~ normal(0, 2);

  // likelihood
  for (j in 1:J) {
    alpha_ln[j] = 0;

    for (p in 1:P_alpha) {
      alpha_ln[j] = alpha_ln[j] + b_alpha[p] * X_alpha[j, p];
    }

    alpha[j] = exp(alpha_ln[j]);

    beta_ln[j] = 0;

    for (p in 1:P_beta) {
      beta_ln[j] = beta_ln[j] + b_beta[p] * X_beta[j, p];
    }

    beta[j] = exp(beta_ln[j]);
  }

  for (i in 1:N) {
    mu[i] = m0[i] - (time[i] / beta[sp[i]])^alpha[sp[i]];
  }

  mT ~ normal(mu, sigma_obs);
}

generated quantities {

  real neg_loglik;
  vector[N_test] mT_pred;

  {
  real alpha_pred;
  real alpha_ln_pred;
  real beta_pred;
  real beta_ln_pred;
  vector[N_test] loglik;

  alpha_ln_pred = 0;

  for (p in 1:P_alpha) {
    alpha_ln_pred = alpha_ln_pred + b_alpha[p] * X_alpha_test[1, p];
  }

  alpha_pred = exp(alpha_ln_pred);

  beta_ln_pred = 0;

  for (p in 1:P_beta) {
    beta_ln_pred = beta_ln_pred + b_beta[p] * X_beta_test[1, p];
  }

  beta_pred = exp(beta_ln_pred);

  for (i in 1:N_test) {
    mT_pred[i] = normal_rng(m0_test[i] - (time_test[i] / beta_pred)^alpha_pred, sigma_obs);
    loglik[i] = normal_lpdf(mT_test[i] | m0_test[i] - (time_test[i] / beta_pred)^alpha_pred, sigma_obs);
  }

  neg_loglik = -1 * sum(loglik);

  }
}

