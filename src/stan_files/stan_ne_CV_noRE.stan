
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
  matrix[J, P] X_k;         // design matrix for alpha effects
  matrix[1, P] X_k_test;
  int<lower=1, upper=J> sp[N];            // species id
}

parameters {
  vector[P] b;
  real<lower=0> sigma_obs;                // observation sd
}

model {

  real k[J];
  real k_ln[J];
  vector[N] k_i;

  // priors
  sigma_obs ~ normal(0, 2);
  b ~ normal(0, 2);

  // likelihood
  for (j in 1:J) {
    k_ln[j] = 0;

    for (p in 1:P) {
      k_ln[j] = k_ln[j] + b[p] * X_k[j, p];
    }

    k[j] = exp(k_ln[j]);
  }

  for (i in 1:N) {
    k_i[i] = k[sp[i]];
  }

  mT ~ normal(m0 - k_i .* time, sigma_obs); // element-wise multiplication of vectors
}

generated quantities {

  real neg_loglik;
  vector[N_test] mT_pred;

  {
  real k_pred;
  real k_ln_pred;
  vector[N_test] loglik;

  k_ln_pred = 0;

  for (p in 1:P) {
    k_ln_pred = k_ln_pred + b[p] * X_k_test[1, p];
  }

  k_pred = exp(k_ln_pred);

  for (i in 1:N_test) {
    mT_pred[i] = normal_rng(m0_test[i] - (k_pred * time_test[i]), sigma_obs);
    loglik[i] = normal_lpdf(mT_test[i] | m0_test[i] - (k_pred * time_test[i]), sigma_obs);
  }

  neg_loglik = -1 * sum(loglik);
  }
}
