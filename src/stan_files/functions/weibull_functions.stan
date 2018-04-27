/**
  * Create weibull model for data
  *
  * @param N number of test data
  * @param m0 initial mass of test data
  * @param time time of test data
  * @param beta predicted beta
  * @param sp list of species' numbers
  * @param alpha predicted alpha
  * @param sigma_obs datapoint error
  * @return mu
  */
vector weibull(int N, vector m0, vector time, vector beta, vector alpha, int[] sp, real sigma_obs) {

  vector[N] mu;

  for (i in 1:N) {
    mu[i] = m0[i] - (time[i] / beta[sp[i]])^alpha[sp[i]];
  }

  return mu;
}

/**
  * Create weibull fit for data
  *
  * @param N number of test data
  * @param m0 initial mass of test data
  * @param time time of test data
  * @param beta_fit fit beta
  * @param alpha_fit fit alpha
  * @param sp species' number
  * @param sigma_obs datapoint error
  * @return A vector of parameter estimates
  */
vector weibull_fit_rng(int N, vector m0, vector time, vector beta_fit, vector alpha_fit, int[] sp, real sigma_obs) {

  vector[N] mT_fit;

  for (i in 1:N) {
    mT_fit[i] = normal_rng(m0[i] - (time[i] / beta_fit[sp[i]])^alpha_fit[sp[i]], sigma_obs);
  }

  return mT_fit;
}

/**
  * Create weibull predictions for test datapoints
  *
  * @param N_test number of test data
  * @param m0_test initial mass of test data
  * @param time_test time of test data
  * @param beta_pred predicted beta
  * @param alpha_pred predicted alpha
  * @param sigma_obs datapoint error
  * @return A vector of parameter estimates
  */
vector weibull_pred_rng(int N_test, vector m0_test, vector time_test, real beta_pred, real alpha_pred, real sigma_obs) {

  vector[N_test] mT_pred;

  for (i in 1:N_test) {
    mT_pred[i] = normal_rng(m0_test[i] - (time_test[i] / beta_pred)^alpha_pred, sigma_obs);
  }

  return mT_pred;
}

/**
  * Calculate negative loglikelihood
  *
  * @param mT_test true mass at time T of test data
  * @param m0_test initial mass of test data
  * @param time_test time of test data
  * @param beta_pred predicted beta
  * @param alpha_pred predicted alpha
  * @param sigma_obs datapoint error
  * @param N_test number of test data
  * @return A vector of parameter estimates
  */
real weibull_negloglik(vector mT_test, vector m0_test, vector time_test, real beta_pred, real alpha_pred, real sigma_obs, int N_test) {

  vector[N_test] loglik;
  real neg_loglik;

  for (i in 1:N_test) {
    loglik[i] = normal_lpdf(mT_test[i] | m0_test[i] - (time_test[i] / beta_pred)^alpha_pred, sigma_obs);
  }

  neg_loglik = -1 * sum(loglik);
  return neg_loglik;
}
