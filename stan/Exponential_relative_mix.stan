// exponential survival mixture cure model
// relative survival

// user-defined functions ----
functions {
  // log hazard
  real log_h (real t, real rate) {
    real logh;
    logh = log(rate);
    return logh;
  }

  real haz (real t, real rate) {
    real h;
    h = rate;
    return h;
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

  // sampling distributions

  real surv_exp_pdf (real t, real d, real rate) {
    real lik;
    lik = haz(t, rate)^d * Surv(t,rate);
    return lik;
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
  // real mu_bg;
  // vector<lower=0> [H] sigma_beta;   // sds of the covariates coefficients
  real<lower=0> sigma_beta;   // sds of the covariates coefficients
  // real<lower=0> sigma_bg;
  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;
  vector[n] h_bg;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  // vector[H] beta_bg;
  real<lower=0, upper=1> curefrac;
}

transformed parameters {
  vector[n] linpred0;
  vector[n] lambda0;
  // vector[n] linpred_bg;
  // vector[n] lambda_bg;

  linpred0 = X*beta0;
  // linpred_bg = X*beta_bg;

  //TODO: why is this not vectorised like everything else?
  for (i in 1:n) {
    lambda0[i] = exp(linpred0[i]);     // rate parameters
    // lambda_bg[i] = exp(linpred_bg[i]);
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);
  // beta_bg ~ normal(mu_bg, sigma_bg);

  curefrac ~ beta(a_cf, b_cf);

  for (i in 1:n) {

    // _known_ point estimate for background survival
    // target += log_mix(curefrac,
    //                   surv_exp_lpdf(t[i] | d[i], h_bg[i]),
    //                   surv_exp_lpdf(t[i] | d[i], h_bg[i] + lambda0[i]));
    // equivalently
    target += log_sum_exp(log(curefrac)
                          + surv_exp_lpdf(t[i] | d[i], h_bg[i]),
                          log1m(curefrac)
                          + surv_exp_lpdf(t[i] | d[i], h_bg[i] + lambda0[i]));


    // // background survival with uncertainty
    // target += log_mix(curefrac,
    //                   surv_exp_lpdf(t[i] | d[i], lambda_bg[i]),
    //                   surv_exp_lpdf(t[i] | d[i], lambda_bg[i] + lambda0[i]));
  }
}

generated quantities {
  real rate0;
  // real rate1;

  rate0 = exp(beta0[1]);
  // rate1 = exp(beta_bg[1]);
}

