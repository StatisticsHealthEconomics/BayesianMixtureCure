// Exponential survival model

functions {
  // Defines the log hazard
  vector log_h (vector t, vector rate) {
    vector[num_elements(t)] logh;
    logh = log(rate);
    return logh;
  }

  // Defines the log survival
  vector log_S (vector t, vector rate) {
    vector[num_elements(t)] logS;
    logS = -rate .* t;
    return logS;
  }

  // Defines the sampling distribution
  real surv_exponential_lpdf (vector t, vector d, vector rate) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,rate) + log_S(t,rate);
    prob = sum(log_lik);
    return prob;
  }
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1 = observed, 0 = censored)
  int H;                  // number of covariates
  matrix[n,H] X;          // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	    // means of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sds of the covariates coefficients

}

parameters {
  vector[H] beta;         // coefficients in the linear predictor (including intercept)
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;

  linpred = X*beta;

  //TODO: why is this not vectorised like everything else?
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);  // rate parameter
  }
}

model {
  beta ~ normal(mu_beta, sigma_beta);
  t ~ surv_exponential(d,mu);
}

generated quantities {
  real rate;                 // rate parameter
  rate = exp(beta[1]);
}

