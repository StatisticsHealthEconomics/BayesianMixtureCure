// generalised gamma survival mixture cure model
// relative survival

functions {

  //  log hazard
  real exp_log_h (real t, real rate) {
    real logh;
    logh = log(rate);
    return logh;
  }

  // exponential distribution hazard
  real exp_haz (real t, real rate) {
    real h;
    h = rate;
    return h;
  }

  // exponential distribution log survival
  // inbuilt exponential_lccdf(y | beta)
  real exp_log_S (real t, real rate) {
    real logS;
    logS = -rate * t;
    return logS;
  }

  real gengamma_Surv (real t, real rate) {
    real Surv;
    real w = (log(t) - mu) / sigma;
    real qq = 1/(Q * Q);                    // shape (gamma)
    real expnu = exp(abs(Q) * w) * qq;      // u

    if (Q == 0) {
      Surv = 1 - normal_cdf(w, 0, 1);
    } else {
      Surv =  1 - gamma_cdf(expnu, qq, 1);
    }
    return Surv;
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

  real surv_gen_gamma_lpdf(real t, real d, real mu, real sigma, real Q) {
    // rescale the distribution accounting for right censoring
    real log_lik;
    real w;
    real tr;
    tr = t * d;
    w = (log(tr) - mu)/sigma;
    log_lik = log(d) - log(sigma*tr) + log(fabs(Q)) + pow(Q, -2)*log(pow(Q, -2)) + pow(Q, -2)*(Q*w - exp(Q*w)) - lgamma(pow(Q, -2));
    return log_lik;
  }

  real joint_exp_gengamma_lpdf(real t, real d, real mu, real scale, real Q, real rate) {
    real log_lik;
    log_lik = d * log(exp_haz(t, rate) + gen_gamma_haz(t, mu, scale, Q)) +
              exp_log_S(t, rate) + gen_gamma_log_S(t, mu, scale, Q);
    return log_lik;
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
  vector[H] mu_bg;
  vector<lower=0> [H] sigma_bg;

  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;
  // vector[n] h_bg;

  real a_Q;                   // generalised gamma hyper-parameters
  real<lower=0> b_Q;
  real a_scale;
  real<lower=0> b_scale;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  vector[H] beta_bg;
  real<lower=0, upper=1> curefrac;
  real Q;
  real<lower=0> scale;
}

transformed parameters {
  vector[n] linpred0;
  vector[n] linpred_bg;
  vector[n] lambda0;
  vector[n] lambda_bg;

  linpred0 = X*beta0;
  linpred_bg = X*beta_bg;

  // rate parameters
  for (i in 1:n) {
    mu[i] = linpred0[i];
    lambda_bg[i] = exp(linpred_bg[i]); // background survival with uncertainty
    // lambda_bg[i] = h_bg[i];           // _known_ point estimate for background survival
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);
  beta_bg ~ normal(mu_bg, sigma_bg);

  scale ~ lognormal(a_scale, b_scale);
  Q ~ normal(a_Q, b_Q);

  curefrac ~ beta(a_cf, b_cf);

  for (i in 1:n) {
    target += log_sum_exp(log(curefrac)
              + surv_exp_lpdf(t[i] | d[i], lambda_bg[i]),
              log1m(curefrac)
              + joint_exp_gengamma_lpdf(t[i] | d[i], mu[i], scale, Q, lambda_bg[i]);
  }
}

generated quantities {
  real m0;
  real rate_bg;
  vector[60] S_bg;
  vector[60] S_0;
  vector[60] S_pred;

  mu0 = beta0[1];
  rate_bg = exp(beta_bg[1]);

  for (i in 1:60) {
    S_bg[i] = exp_Surv(i, rate_bg);
    S_0[i] = gengamma_Surv(i, mu0, sigma, Q);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S_bg[i]*S_0[i];
  }
}

