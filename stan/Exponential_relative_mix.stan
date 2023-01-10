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

  real exp_pdf (real t, real rate) {
    real lik;
    lik = rate*exp(-rate*t);
    return lik;
  }

  real exp_lpdf (real t, real rate) {
    real log_lik;
    log_lik = log(rate) - rate*t;
    return log_lik;
  }

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

  // intercept only -
  // real mu_beta;	              // means of the covariates coefficients
  // real mu_bg;
  // real<lower=0> sigma_beta;    // sds of the covariates coefficients
  // real<lower=0> sigma_bg;
  // intercept and gradient -
  vector[H] mu_beta;
  vector<lower=0> [H] sigma_beta;

  // vector[H] mu_bg;
  // vector<lower=0> [H] sigma_bg;
  vector[n] h_bg;

  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  // vector[H] beta_bg;
  real<lower=0, upper=1> curefrac;  //TODO: define as simplex?
}

transformed parameters {
  vector[n] linpred0;
  // vector[n] linpred_bg;
  vector[n] lambda0;
  vector[n] lambda_bg;

  linpred0 = X*beta0;
  // linpred_bg = X*beta_bg;

  //TODO: is there a vectorised exp?
  // rate parameters
  for (i in 1:n) {
    lambda0[i] = exp(linpred0[i]);
    // lambda_bg[i] = exp(linpred_bg[i]); // background survival with uncertainty
    lambda_bg[i] = h_bg[i];           // _known_ point estimate for background survival
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);
  // beta_bg ~ normal(mu_bg, sigma_bg);

  curefrac ~ beta(a_cf, b_cf);

  for (i in 1:n) {

    // joint survival
    target += log_sum_exp(log(curefrac)
                        + log_S(t[i], lambda_bg[i]),
                        log1m(curefrac)
                        + log_S(t[i], lambda_bg[i]) + log_S(t[i], lambda0[i]));
    // joint hazard
    target += d[i] * log_sum_exp(log(lambda_bg[i]),
                                 log(1 - curefrac) + exp_lpdf(t[i], lambda0[i]) -
                                  log(curefrac + (1 - curefrac)*Surv(t[i], lambda0[i])));
  }
}

generated quantities {
  real rate0;
  real rate_bg;
  vector[60] S_bg;
  vector[60] S_0;
  vector[60] S_pred;

  rate0 = exp(beta0[1]);
  // rate_bg = exp(beta_bg[1]);
  rate_bg = mean(h_bg);

  for (i in 1:60) {
    S_bg[i] = Surv(i, rate_bg);
    S_0[i] = Surv(i, rate_bg + rate0);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S_0[i];
  }
}

