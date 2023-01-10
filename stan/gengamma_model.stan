// generalised gamma survival  model

functions {

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

##TODO: check with survHE is this right?
## https://github.com/n8thangreen/survHE/blob/main/src/stan_files/GenGamma.stan
##      is censoring correct? is there are log(0)??
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

## this should give the same results as surv_gen_gamma_lpdf. does it?
  real gengamma_lpdf(real t, real d, real mu, real scale, real Q) {
    real log_lik;
    log_lik = d * log(gen_gamma_haz(t, mu, scale, Q)) + gen_gamma_log_S(t, mu, scale, Q);
    return log_lik;
  }
}

data {
  int<lower=0> n;             // number of observations
  vector[n] t;                // observed times
  int H;                      // number of covariates
  matrix[n,H] X;              // matrix of covariates (with n rows and H columns)

  vector[H] mu_beta;
  vector<lower=0> [H] sigma_beta;

  real a_Q;                   // generalised gamma hyper-parameters
  real<lower=0> b_Q;
  real a_scale;
  real<lower=0> b_scale;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  real Q;
  real<lower=0> scale;
}

transformed parameters {
  vector[n] linpred0;
  vector[n] mu;

  linpred0 = X*beta0;

  for (i in 1:n) {
    mu[i] = linpred0[i];
  }
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);

  scale ~ lognormal(a_scale, b_scale);
  Q ~ normal(a_Q, b_Q);

  for (i in 1:n) {
    target += gen_gamma_lpdf(t[i] | mu[i], scale, Q);
  }
}

generated quantities {
  real mu0;
  real rate_bg;
  vector[60] S_pred;

  mu0 = beta0[1];

  for (i in 1:60) {
    S_pred[i] = gen_gamma_Surv(i, mu0, scale, Q);
  }
}

