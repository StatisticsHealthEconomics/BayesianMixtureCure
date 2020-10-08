
# run stan mixture cure model
# real-life dataset
#

library(rstan)
library(shinystan)
library(dplyr)

load("C:/Users/Nathan/Documents/R/mixture_cure_model/data/bg_hazards.RData")

tx_dat <-
  bg_hazards %>%
  select(TRTA, pfs, pfs_event, PFS_rate) %>%
  mutate(PFS_rate = ifelse(PFS_rate == 0, 0.00001, PFS_rate)) %>%
  split(bg_hazards$TRTA)

tx_name <- "IPILIMUMAB"

data_list <-
  list(
    n = nrow(tx_dat[[tx_name]]),
    t = tx_dat[[tx_name]]$pfs,
    d = tx_dat[[tx_name]]$pfs_event,
    H = 1,
    X = matrix(rep(1, nrow(tx_dat[[tx_name]])), ncol = 1),
    mu_beta = 0,
    sigma_beta = 1,
    a_cf = 1,
    b_cf = 1,
    h_bg = tx_dat[[tx_name]]$PFS_rate)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# stan_rdump(c("n_obs", "y"), file = "mix.data.R")

n_iter <- 2000#0
n_warmup <- 1000#0

stan_base <-
  stan(
    file = here::here("stan", "Exponential_relative_mix.stan"),
    data = data_list,
    warmup = n_warmup,
    # control = list(adapt_delta = 0.9,
    #                max_treedepth = 20),
    iter = n_iter,
    chains = 1)

stan_base

