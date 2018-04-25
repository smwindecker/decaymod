
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

  vector[N] mT_fit;

  {
    real alpha_fit[J];
    real beta_fit[J];
    real alpha_ln_fit[J];
    real beta_ln_fit[J];

    for (j in 1:J) {
      alpha_ln_fit[j] = 0;

      for (p in 1:P_alpha) {
        alpha_ln_fit[j] = alpha_ln_fit[j] + b_alpha[p] * X_alpha[j, p];
      }

      alpha_fit[j] = exp(alpha_ln_fit[j]);

      beta_ln_fit[j] = 0;

      for (p in 1:P_beta) {
        beta_ln_fit[j] = beta_ln_fit[j] + b_beta[p] * X_beta[j, p];
      }

      beta_fit[j] = exp(beta_ln_fit[j]);
    }

    for (i in 1:N) {
      mT_fit[i] = normal_rng(m0[i] - (time[i] / beta_fit[sp[i]])^alpha_fit[sp[i]], sigma_obs);
    }

  }
}

