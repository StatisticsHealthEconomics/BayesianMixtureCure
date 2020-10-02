
#
# script for running stan
# survival models


library(rstan)
library(shinystan)

# data <- read.csv(here::here("raw_data", "data.csv"), header = TRUE)

n_obs <- 10
times0 <- rexp(n_obs, 2)
times1 <- rexp(n_obs, 3)

data_list <-
  list(
    n = 2*n_obs,                 # number of observations
    t = c(times0, times1),       # observed times
    d = rep(1, 2*n_obs),         # censoring indicator (1 = observed, 0 = censored)
    H = 1,                       # number of covariates
    X = matrix(c(rep(0, n_obs),  # matrix of covariates (with n rows and H columns)
                 rep(1, n_obs)),
               ncol = 1),
    mu_beta = 0,	               # mean of the covariates coefficients
    sigma_beta = 1               # sd of the covariates coefficients
  )


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter <- 200#000
n_warmup <- 100#000

stan_base <-
  stan(
    file = here::here("stan", "Exponential.stan"),
    data = data_list,
    warmup = n_warmup,
    # control = list(adapt_delta = 0.9,
    #                max_treedepth = 20),
    iter = n_iter,
    chains = 1)

# launch_shinystan(stan_base)
# base <- extract(stan_base)

