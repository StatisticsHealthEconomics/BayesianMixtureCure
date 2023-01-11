// generalised gamma survival mixture cure model
// splitting up the censored and non-censored likelihood parts

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

  real exp_pdf(real t, real rate) {
    real x;
    x = rate*exp(-rate * t);
    return x;
  }

  real exp_lpdf(real t, real rate) {
    real x;
    x = log(rate) - (rate * t);
    return x;
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
  // generalised gamma

  real gen_gamma_Surv (real t, real mu, real sigma, real Q) {
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

  real gen_gamma_log_S(real t, real mu, real sigma, real Q) {
    real log_S;
    log_S = log(gen_gamma_Surv(t, mu, sigma, Q));
    return log_S;
  }

  real gen_gamma_lpdf(real t, real mu, real sigma, real Q) {
    real log_lik;
    real w;
    w = (log(t) - mu)/sigma;
    log_lik = -log(sigma*t) + log(fabs(Q)) + pow(Q, -2)*log(pow(Q, -2)) + pow(Q, -2)*(Q*w-exp(Q*w)) - lgamma(pow(Q, -2));
    return log_lik;
  }

  real gen_gamma_pdf(real t, real mu, real sigma, real Q) {
    real x;
    x = exp(gen_gamma_lpdf(t | mu, sigma, Q));
    return x;
  }

  real gen_gamma_log_h(real t, real mu, real sigma, real Q) {
    real log_h;
    log_h = gen_gamma_lpdf(t | mu, sigma, Q) - gen_gamma_log_S(t, mu, sigma, Q);
    return log_h;
  }

  real gen_gamma_haz(real t, real mu, real sigma, real Q) {
    real haz;
    haz = exp(gen_gamma_log_h(t, mu, sigma, Q));
    return haz;
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
  int<lower=0> N;             // number of observations
  int<lower=0> n_cens;
  int<lower=0> n_unc;
  vector[n_cens] t_cens;           // observed times
  vector[n_unc] t_unc;
  int<lower=1> cens_idx[n_cens];      // censored index
  int<lower=1> unc_idx[n_unc];
  int H;                      // number of covariates
  matrix[N,H] X;              // matrix of covariates (with n rows and H columns)

  vector[H] mu_beta;
  vector<lower=0> [H] sigma_beta;

  vector[N] h_bg;             // fixed hazard

  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;

  real a_Q;                   // generalised gamma hyper-parameters
  real<lower=0> b_Q;
  real a_scale;
  real<lower=0> b_scale;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  real<lower=0, upper=1> curefrac;
  real Q;
  real<lower=0> scale;
}

transformed parameters {
  vector[N] linpred0;
  vector[N] mu;
  vector[N] lambda_bg;

  linpred0 = X*beta0;

  for (i in 1:N) {
    mu[i] = linpred0[i];
    lambda_bg[i] = h_bg[i];           // _known_ point estimate for background survival
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);

  scale ~ lognormal(a_scale, b_scale);
  Q ~ normal(a_Q, b_Q);

  curefrac ~ beta(a_cf, b_cf);

  if (n_cens > 0) {
    for (i in 1:n_cens) {
      target += log_sum_exp(log(curefrac)
                        + exp_log_S(t_cens[i], lambda_bg[cens_idx[i]]),
                        log1m(curefrac)
                        + exp_log_S(t_cens[i], lambda_bg[cens_idx[i]]) + gen_gamma_log_S(t_cens[i], mu[cens_idx[i]], scale, Q));
    }
  }

  for (i in 1:n_unc) {
    target += log_sum_exp(log(curefrac)
                        + exp_lpdf(t_unc[i] | lambda_bg[unc_idx[i]]),
                        log1m(curefrac)
                        + log(gen_gamma_pdf(t_unc[i], mu[unc_idx[i]], scale, Q)*exp_Surv(t_unc[i], lambda_bg[unc_idx[i]]) +
                              exp_pdf(t_unc[i], lambda_bg[unc_idx[i]])*gen_gamma_Surv(t_unc[i], mu[unc_idx[i]], scale, Q)));
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
    S_0[i] = gen_gamma_Surv(i, mu0, scale, Q);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S_bg[i]*S_0[i];
  }
}

