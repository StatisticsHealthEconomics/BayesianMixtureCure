
# run stan mixture cure model
# with simulated data
# generalised gamma distribution


library(rstan)
library(shinystan)
library(dplyr)
library(flexsurv)
library(survival)

n <- 1000
x <- rnorm(n)
loc <- 2 + 0.2*x
cf <- 2/14

dat <- data.frame(
  t_gg = flexsurv::rgengamma(n, mu = loc, sigma = log(1.5), Q = -0.2),
  t_exp = rexp(n, rate = 0.001),
  cure = runif(n) < cf) |>
  mutate(t = ifelse(cure, t_exp, pmin(t_exp, t_gg)),
         d = 1,
         x = as.numeric(x))

survfit(Surv(t, d)~1, data = dat) |> plot(xlim = c(0,50))
survfit(Surv(t_gg, d)~1, data = dat) |> lines(col = "red")
survfit(Surv(t_exp, d)~1, data = dat) |> lines(col = "green")

data_list <-
  list(
    n = n,
    t = dat$t,
    d = dat$d,
    H = 2,
    X = matrix(c(rep(1, n), x),
               byrow = FALSE,
               ncol = 2),
    mu_beta = c(0.3, 0),
    sigma_beta = c(0.8, 0.2),
    a_Q = -0.2,
    b_Q = 0.6,
    a_scale = log(1.5),
    b_scale = 0.1,
    a_cf = 2,
    b_cf = 12,
    h_bg = rep(0.001, n))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter <- 500#0
n_warmup <- 100#0

stan_base <-
  stan(
    file = here::here("stan", "gengamma_mixture_cure_model.stan"),
    data = data_list,
    control = list(adapt_delta=0.99, max_treedepth=11),
    warmup = n_warmup,
    iter = n_iter,
    chains = 1)

stan_base

