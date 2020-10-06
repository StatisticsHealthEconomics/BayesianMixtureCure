// exponential survival mixture cure model
//
// see https://mc-stan.org/docs/2_19/stan-users-guide/summing-out-the-responsibility-parameter.html


// user-defined functions ----
functions {
  // log hazard
  real log_h (real t, real rate) {
    real logh;
    logh = log(rate);
    return logh;
  }

  // log survival
  real log_S (real t, real rate) {
    real logS;
    logS = -rate * t;
    return logS;
  }

  real Surv (real t, real rate) {
    real S;
    S = exp(-rate * t);
    return S;
  }

  real haz (real t, real rate) {
    real h;
    h = rate;
    return h;
  }

  // sampling distributions

  real surv_exp_pdf (real t, real d, real rate) {
    real lik;
    lik = haz(t, rate)^d + Surv(t,rate);
    return lik;
  }

  //TODO:
  real surv_exp_mix_lpdf (vector t, vector d, vector rate0, vector rate1, real curefrac) {
    vector[num_elements(t)] log_lik;
    real prob;
    // log_lik = log(exp(log(curefrac) + surv_exp_log_lik(t, d, rate0)) + exp(log(1 - curefrac) + surv_exp_log_lik(t, d, rate1)));
    for (i in 1:num_elements(t)) {
      log_lik[i] = log(curefrac * surv_exp_pdf(t[i], d[i], rate0[i]) + (1 - curefrac) * surv_exp_pdf(t[i], d[i], rate1[i]));
    }
    prob = sum(log_lik);
    return prob;
  }

  real surv_exp_lpdf (real t, real d, real rate) {
    real log_lik;
    log_lik = d * log_h(t, rate) + log_S(t, rate);
    return log_lik;
  }
}

// input data ----
data {
  int<lower=0> n;             // number of observations
  vector[n] t;                // observed times
  vector[n] d;                // censoring indicator (1 = observed, 0 = censored)
  int H;                      // number of covariates
  matrix[n,H] X;              // matrix of covariates (with n rows and H columns)
  // vector[H] mu_beta;	      // means of the covariates coefficients
  real mu_beta;	              // means of the covariates coefficients
  // vector<lower=0> [H] sigma_beta;   // sds of the covariates coefficients
  real<lower=0> sigma_beta;   // sds of the covariates coefficients
  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  vector[H] beta1;
  real<lower=0, upper=1> curefrac;
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
    mu1[i] = exp(linpred1[i]);
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);
  beta1 ~ normal(mu_beta, sigma_beta);
  curefrac ~ beta(a_cf, b_cf);

  for (i in 1:n)
    target += log_mix(curefrac,
                      surv_exp_lpdf(t[i] | d[i], mu0[i]),
                      surv_exp_lpdf(t[i] | d[i], mu1[i]));
                      // surv_exp_lpdf(t[i] | d[i], mu0[i]) + surv_exp_lpdf(t[i] | d[i], mu1[i]));

  // target += log_sum_exp(log(curefrac) + surv_exp_lpdf(t[i] | d, mu0),
  //                       log(1 - curfrac) + surv_exp_lpdf(t[i] | d, mu0) + surv_exp_lpdf(t | d, mu1));
}

generated quantities {
  real rate0;
  real rate1;
  rate0 = exp(beta0[1]);
  rate1 = exp(beta1[1]);
}

