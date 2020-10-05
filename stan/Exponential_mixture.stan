// Exponential survival model

functions {
  // define the log hazard
  vector log_h (vector t, vector rate) {
    vector[num_elements(t)] logh;
    logh = log(rate);
    return logh;
  }

  // define the log survival
  vector log_S (vector t, vector rate) {
    vector[num_elements(t)] logS;
    logS = -rate .* t;
    return logS;
  }

  // define the sampling distribution
  real surv_exponential_lpdf (vector t, vector d, vector rate0, vector rate1, vector curefrac) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = curefrac .*(d .* log_h(t,rate0) + log_S(t,rate0)) + (1 - curefrac) .*(d .* log_h(t,rate1) + log_S(t,rate1));
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
  vector[1] a_cf;
  vector[1] b_cf;
}

parameters {
  vector[H] beta0;         // coefficients in the linear predictor (including intercept)
  vector[H] beta1;
  vector[1] curefrac;
}

transformed parameters {
  vector[n] linpred0;
  vector[n] mu0;
  vector[n] linpred1;
  vector[n] mu1;

  linpred0 = X*beta0;
  linpred1 = X*beta1;

  //TODO: why is this not vectorised like everything else?
  for (i in 1:n) {
    mu0[i] = exp(linpred0[i]);  // rate parameter
    mu1[i] = exp(linpred1[i]);  // rate parameter
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);
  beta1 ~ normal(mu_beta, sigma_beta);
  curefrac ~ beta(a_cf, b_cf);
  t ~ surv_exponential(d, mu0, mu1, curefrac);
}

generated quantities {
  real rate0;
  real rate1;
  rate0 = exp(beta0[1]);
  rate1 = exp(beta1[1]);
}

