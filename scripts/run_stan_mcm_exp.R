
#
# script for running stan
# survival models

##use model.amtrix()?

library(rstan)
library(shinystan)

# data <- read.csv(here::here("raw_data", "data.csv"), header = TRUE)

n_obs <- 20
p_cure <- 0.2
n_cure <- n_obs*p_cure

times_c <- rexp(n_cure, 0.5)
times_n <- rexp(n_obs - n_cure, 2)

data_list <-
  list(
    n = n_obs,                   # number of observations
    t = c(times_c, times_n),     # observed times
    d = rep(1, n_obs),           # censoring indicator (1 = observed, 0 = censored)
    H = 1,                       # number of covariates
    X = rep(1, n_obs),           # matrix of covariates (with n rows and H columns)
    mu_beta = 0,	               # mean of the covariates coefficients
    sigma_beta = 1,              # sd of the covariates coefficients
    a_cf = 1,                    # cure fraction beta distn
    b_cf = 1
  )

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter <- 200000
n_warmup <- 100000

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


## plots

library(survHE)

model_formula <- as.formula(Surv(t, d) ~ group)

data_df <- data.frame(d = data_list$d,
                      t = data_list$t,
                      group = group)

my_hmc <-
  list(models = list(Exponential = stan_base),
       model.fitting = NULL,
       method = "hmc",
       misc = list(data = data_df,
                   vars = list(time = "t",
                               event = "d",
                               factors = "group",
                               covs = NULL,
                               nlevs = 2),
                   formula = model_formula,
                   km = rms::npsurv(formula = model_formula,
                                    data = data_df),
                   data.stan = list(data_list))
  )

class(my_hmc) <- "survHE"

# hmc <-
#   fit.models(
#     formula = Surv(recyrs, censrec) ~ group,
#     data = bc,
#     distr = "exp",
#     method = "hmc")
# hmc_surv <- make.surv(hmc)
# plot(hmc)

my_surv <- make.surv(my_hmc)

plot(my_hmc)

# save()

# data <- data.frame(x=rexp(n = 100000, rate = 2))
# m <- ggplot(data, aes(x=data$x))
# m + geom_density() + xlim(0,10)


