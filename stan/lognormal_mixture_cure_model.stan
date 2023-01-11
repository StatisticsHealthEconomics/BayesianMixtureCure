// generalised gamma survival mixture cure model
// relative survival

functions {

  ////////////////
    // exponential

  real exp_log_h (real t, real rate) {
    real logh;
    logh = log(rate);
    return logh;
  }

  real exp_haz (real t, real rate) {
    real h;
    h = rate;
    return h;
  }

  // inbuilt exponential_lccdf(y | beta)
  real exp_log_S (real t, real rate) {
    real logS;
    logS = -rate * t;
    return logS;
  }

  real exp_Surv(real t, real rate) {
    real S;
    S = exp(-rate * t);
    return S;
  }

  real surv_exp_lpdf (real t, real d, real rate) {
    real log_lik;
    log_lik = d * exp_log_h(t, rate) + exp_log_S(t, rate);
    return log_lik;
  }

  //////////////////////
  // log normal

  real lognormal_log_S (real t, real mu, real sigma) {
    real log_S;
    log_S = log(1 - Phi((log(t) - mu)/sigma));
    return log_S;
  }

  real lognormal_Surv (real t, real mu, real sigma) {
   real Surv;
   Surv = 1 - Phi((log(t) - mu)/sigma);
   return Surv;
  }
}

data {
  int<lower=0> n;             // number of observations
  vector[n] t;                // observed times
  vector[n] d;                // censoring indicator (1 = observed, 0 = censored)
  int H;                      // number of covariates
  matrix[n,H] X;              // matrix of covariates (with n rows and H columns)

  vector[H] mu_beta;
  vector<lower=0> [H] sigma_beta;

  vector[n] h_bg;             // fixed hazard

  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;

  real a_scale;
  real<lower=0> b_scale;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  real<lower=0, upper=1> curefrac;
  real<lower=0> scale;
}

transformed parameters {
  vector[n] linpred0;
  vector[n] mu;
  vector[n] lambda_bg;

  linpred0 = X*beta0;

  for (i in 1:n) {
    mu[i] = linpred0[i];
    lambda_bg[i] = h_bg[i];           // _known_ point estimate for background survival
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);

  scale ~ lognormal(a_scale, b_scale);

  curefrac ~ beta(a_cf, b_cf);

  for (i in 1:n) {
    target += log_sum_exp(log(curefrac)
                          + exp_log_S(t[i], lambda_bg[i]),
                          log1m(curefrac)
                          + exp_log_S(t[i], lambda_bg[i]) + lognormal_log_S(t[i], mu[i], scale));

    target += d[i] * log_sum_exp(log(lambda_bg[i]),
                                 log(1 - curefrac) + lognormal_lpdf(t[i] | mu[i], scale) -
                                   log(curefrac + (1 - curefrac)*lognormal_Surv(t[i], mu[i], scale)));

// all stan core functions
    // target += log_sum_exp(log(curefrac)
    //                       + exponential_lccdf(t[i] | lambda_bg[i]),
    //                       log1m(curefrac)
    //                       + exponential_lccdf(t[i] | lambda_bg[i]) + lognormal_lccdf(t[i] | mu[i], scale));
    //
    // target += d[i] * log_sum_exp(log(lambda_bg[i]),
    //                              log(1 - curefrac) + lognormal_lpdf(t[i] | mu[i], scale) -
    //                                log(curefrac + (1 - curefrac)*(1-lognormal_cdf(t[i] | mu[i], scale)))));
  }
}

generated quantities {
  real mu0;
  real rate_bg;
  vector[60] S_bg;
  vector[60] S_0;
  vector[60] S_pred;

  mu0 = beta0[1];

  rate_bg = mean(h_bg);

  for (i in 1:60) {
    S_bg[i] = exp_Surv(i, rate_bg);
    S_0[i] = lognormal_Surv(i, mu0, scale);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S_bg[i]*S_0[i];
  }
}

