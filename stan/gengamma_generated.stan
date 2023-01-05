// generalised gamma survival mixture cure model
// relative survival

functions {

  ////////////////
  // exponential

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
  real a_cf;
  real b_cf;
  real a_mu;                   // generalised gamma hyper-parameters
  real<lower=0> b_mu;
  real a_Q;
  real<lower=0> b_Q;
  real a_scale;
  real<lower=0> b_scale;

  // int N_samples;
}

parameters {
}

model {
}

generated quantities {

  real mu;
  real scale;
  real Q;
  real<lower=0, upper=1> cf;

  // real cf = 0.3; //0.1;
  real rate_bg = 0.002;

  vector[60] S_bg;
  vector[60] S_0;
  vector[60] S_pred;

  for (i in 1:60) {
    mu = normal_rng(a_mu, b_mu);
    scale = lognormal_rng(a_scale, b_scale);
    Q = normal_rng(a_Q, b_Q);
    cf = beta_rng(a_cf, b_cf);

    S_bg[i] = exp_Surv(i, rate_bg);
    S_0[i] = gen_gamma_Surv(i, mu, scale, Q);
    S_pred[i] = cf*S_bg[i] + (1 - cf)*S_bg[i]*S_0[i];
  }
}

