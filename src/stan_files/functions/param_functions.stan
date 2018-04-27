/**
  * Create likelihoods for parameter with random effects
  *
  * @param b effects vector of parameters
  * @param X model matrix of effect values
  * @param P number of effects
  * @param J number of species
  * @param a random effect
  * @return A vector of parameter estimates
  */
  vector param_re(vector b, matrix X, int P, int J, vector a) {

    vector[J] param_ln;
    vector[J] param_vector;

    for (j in 1:J) {

      param_ln[j] = a[j];

      for (p in 1:P) {
        param_ln[j] = param_ln[j] + b[p] * X[j, p];
      }

      param_vector[j] = exp(param_ln[j]);
    }

    return param_vector;
  }

/**
  * Create likelihoods for parameter without random effects
  *
  * @param b effects vector of parameters
  * @param X model matrix of effect values
  * @param P number of effects
  * @param J number of species
  * @return A vector of parameter estimates
  */
  vector param(vector b, matrix X, int P, int J) {

    vector[J] param_ln;
    vector[J] param_vector;

    for (j in 1:J) {

      param_ln[j] = 0;

      for (p in 1:P) {
        param_ln[j] = param_ln[j] + b[p] * X[j, p];
      }

      param_vector[j] = exp(param_ln[j]);
    }

    return param_vector;
  }

/**
  * Create prediction for parameter
  *
  * @param b effects vector of parameters
  * @param X model matrix of effect values
  * @param P number of effects
  * @return A vector of parameter estimates
  */
  real param_pred(vector b, matrix X, int P) {

    real param_ln;
    real param_value;

    param_ln = 0;

    for (p in 1:P) {
      param_ln = param_ln + b[p] * X[1, p];
    }

    param_value = exp(param_ln);
    return param_value;
  }
