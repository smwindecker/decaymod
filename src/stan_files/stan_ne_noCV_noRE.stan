
data {
  int<lower=1> N;                         // number of data points
  vector[N] mT;                             // logged mass at time T
  vector[N] m0;                             // logged initial mass
  vector<lower=0>[N] time;                  // time
  int<lower=1> P;                   // no fixefs
  int<lower=1> J;                         // number of species
  matrix[J, P] X_k;         // design matrix for alpha effects
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

  vector[N] mT_fit;

  {
  real k_ln_fit[J];
  real k_fit[J];

  for (j in 1:J) {
    k_ln_fit[j] = 0;

    for (p in 1:P) {
      k_ln_fit[j] = k_ln_fit[j] + b[p] * X_k[j, p];
    }

    k_fit[j] = exp(k_ln_fit[j]);
  }

  for (i in 1:N) {
    mT_fit[i] = normal_rng(m0[i] - (k_fit[sp[i]] * time[i]), sigma_obs);
  }

  }
}
