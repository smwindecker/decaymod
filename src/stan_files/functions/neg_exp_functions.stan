/**
  * Create negative exponential fit for data
  *
  * @param N number of data
  * @param m0 initial mass of data
  * @param time time of data
  * @param k_fit fit k
  * @param sp species' number
  * @param sigma_obs datapoint error
  * @return A vector of parameter estimates
  */
vector negexp_fit_rng(int N, vector m0, vector time, vector k_fit, int[] sp, real sigma_obs) {

  vector[N] mT_fit;

  for (i in 1:N) {
    mT_fit[i] = normal_rng(m0[i] - (k_fit[sp[i]] * time[i]), sigma_obs);
  }

  return mT_fit;
}

/**
  * Create negative exponential fit for simulated data
  *
  * @param N_sim number of simulated data
  * @param m0_sim simulated m0
  * @param time_sim simulated time data
  * @param k_fit fit k
  * @param sp_sim species' number
  * @return A vector of parameter estimates
  */
vector negexp_sim_rng(int N_sim, vector m0_sim, vector time_sim, vector k_fit, int[] sp_sim) {

  vector[N_sim] mT_sim;

  for (i in 1:N_sim) {
    mT_sim[i] = m0_sim[i] - (k_fit[sp_sim[i]] * time_sim[i]);
  }

  return mT_sim;
}

/**
  * Create negative exponential predictions for test datapoints
  *
  * @param N_test number of test data
  * @param m0_test initial mass of test data
  * @param time_test time of test data
  * @param k_pred predicted k
  * @param sigma_obs datapoint error
  * @return A vector of parameter estimates
  */
vector negexp_pred_rng(int N_test, vector m0_test, vector time_test, real k_pred, real sigma_obs) {

  vector[N_test] mT_pred;

  for (i in 1:N_test) {
    mT_pred[i] = normal_rng(m0_test[i] - (k_pred * time_test[i]), sigma_obs);
  }

  return mT_pred;
}

/**
  * Calculate negative loglikelihood for negative exponential
  *
  * @param mT_test true mass at time T of test data
  * @param m0_test initial mass of test data
  * @param time_test time of test data
  * @param k_pred predicted k
  * @param sigma_obs datapoint error
  * @param N_test number of test data
  * @return A vector of parameter estimates
  */
real negexp_negloglik(vector mT_test, vector m0_test, vector time_test, real k_pred, real sigma_obs, int N_test) {

  vector[N_test] loglik;
  real neg_loglik;

  for (i in 1:N_test) {
    loglik[i] = normal_lpdf(mT_test[i] | m0_test[i] - (k_pred * time_test[i]), sigma_obs);
  }

  neg_loglik = -1 * sum(loglik);
  return neg_loglik;
}

