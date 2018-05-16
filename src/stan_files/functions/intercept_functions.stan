
/**
  * Create intercept fit for data
  *
  * @param N number of test data
  * @param m0 initial mass of test data
  * @param time time of test data
  * @param i_fit fit intercept
  * @param k_fit fit k
  * @param sp species' number
  * @param sigma_obs datapoint error
  * @return A vector of parameter estimates
  */
vector intercept_fit_rng(int N, vector m0, vector time, vector i_fit, vector k_fit, int[] sp, real sigma_obs) {

  vector[N] mT_fit;

  for (i in 1:N) {
    mT_fit[i] = normal_rng(m0[i] - (k_fit[sp[i]] * time[i]) + i_fit[sp[i]], sigma_obs);
  }

  return mT_fit;
}

/**
  * Create intercept predictions for test datapoints
  *
  * @param N_test number of test data
  * @param m0_test initial mass of test data
  * @param time_test time of test data
  * @param i_pred predicted intercept
  * @param k_pred predicted k
  * @param sigma_obs datapoint error
  * @return A vector of parameter estimates
  */
vector intercept_pred_rng(int N_test, vector m0_test, vector time_test, real i_pred, real k_pred, real sigma_obs) {

  vector[N_test] mT_pred;

  for (i in 1:N_test) {
    mT_pred[i] = normal_rng(m0_test[i] - (k_pred * time_test[i]) + i_pred, sigma_obs);
  }

  return mT_pred;
}

/**
  * Calculate negative loglikelihood
  *
  * @param mT_test true mass at time T of test data
  * @param m0_test initial mass of test data
  * @param time_test time of test data
  * @param i_pred predicted intercept
  * @param k_pred predicted k
  * @param sigma_obs datapoint error
  * @param N_test number of test data
  * @return A vector of parameter estimates
  */
real intercept_negloglik(vector mT_test, vector m0_test, vector time_test, real i_pred, real k_pred, real sigma_obs, int N_test) {

  vector[N_test] loglik;
  real neg_loglik;

  for (i in 1:N_test) {
    loglik[i] = normal_lpdf(mT_test[i] | m0_test[i] - (k_pred * time_test[i]) + i_pred, sigma_obs);
  }

  neg_loglik = -1 * sum(loglik);
  return neg_loglik;
}

