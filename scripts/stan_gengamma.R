
# run stan model
# generalised gamma distribution


library(rstan)
library(shinystan)
library(dplyr)
library(flexsurv)
library(survival)

n <- 1000
x <- rnorm(n)
loc <- 2 + 0.2*x
dat <- data.frame(
  t = flexsurv::rgengamma(n, mu = loc, sigma = log(1.5), Q = -0.2),
  d = 1,
  x = as.numeric(x))

survfit(Surv(t, d)~1, data=dat) |> plot()
coxph(Surv(t, d)~1+x, data=dat)


data_list <-
  list(
    n = n,
    t = dat$t,
    H = 2,
    X = matrix(c(rep(1, n),
                 dat$x),
               byrow = FALSE,
               ncol = 2),
    mu_beta = c(0.3, 0),
    sigma_beta = c(0.8, 0.2),
    a_Q = -0.2,
    b_Q = 0.6,
    a_scale = log(1.5),
    b_scale = 0.1)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter <- 500#0
n_warmup <- 100#0

stan_base <-
  stan(
    file = here::here("stan", "gengamma_model.stan"),
    data = data_list,
    warmup = n_warmup,
    iter = n_iter,
    chains = 1)

stan_base


########
# plots

library(ggplot2)

fit_stan <- extract(stan_base)

plot_dat <-
  data.frame(s_pred = colMeans(fit_stan$S_pred),
             time = 1:60)

gg <-
  ggplot(plot_dat, aes(x = time, y = s_pred)) +
  geom_line()


km <- survfit(Surv(dat$t, dat$d) ~ 1)
km_data <- data.frame(surv = km$surv,
                      time = km$time)

gg + geom_step(aes(x = time, y = surv),
               linewidth = 1,
               data = km_data,
               inherit.aes = FALSE)

