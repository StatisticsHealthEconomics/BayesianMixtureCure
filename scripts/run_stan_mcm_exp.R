
# run_stan_mcm_exp.R
# script for running stan
# survival mixture cure models


library(rstan)
library(shinystan)


n_obs <- 1000
p_cure <- 0.45
rates <- c(0.5, 2)

z <- rbinom(n_obs, 1, p_cure) + 1
times <- rexp(n_obs, rates[z])

data_list <-
  list(
    n = n_obs,                               # number of observations
    t = times,                               # observed times
    d = rep(1, n_obs),                       # censoring indicator (1 = observed, 0 = censored)
    H = 1,                                   # number of covariates
    X = matrix(rep(1, n_obs), ncol = 1),     # matrix of covariates (with n rows and H columns)
    mu_beta = 0,	                           # mean of the covariates coefficients
    sigma_beta = 1,                          # sd of the covariates coefficients
    a_cf = 1,                                # cure fraction beta distn
    b_cf = 1
  )

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# stan_rdump(c("n_obs", "y"), file = "mix.data.R")

n_iter <- 2000#0
n_warmup <- 1000#0

stan_base <-
  stan(
    file = here::here("stan", "Exponential_mixture.stan"),
    data = data_list,
    warmup = n_warmup,
    # control = list(adapt_delta = 0.9,
    #                max_treedepth = 20),
    iter = n_iter,
    chains = 1)

stan_base

# launch_shinystan(stan_base)
# base <- extract(stan_base)

