
# run stan mixture cure model
# real-life dataset
#

library(rstan)
library(shinystan)
library(dplyr)

load("C:/Users/n8tha/Documents/R/bgfscure/data/surv_input_data.RData")

tx_dat <-
  surv_input_data %>%
  select(TRTA, pfs, pfs_event, PFS_rate, PFSage) %>%
  mutate(PFS_rate =
           ifelse(PFS_rate == 0, 0.00001, PFS_rate)) %>% # replace 0
  split(surv_input_data$TRTA)

tx_name <- "IPILIMUMAB"

data_list <-
  list(
    n = nrow(tx_dat[[tx_name]]),
    t = tx_dat[[tx_name]]$pfs,
    d = tx_dat[[tx_name]]$pfs_event,
    H = 2,
    X = matrix(c(rep(1, nrow(tx_dat[[tx_name]])),
               tx_dat[[tx_name]]$PFSage),
               byrow = FALSE,
               ncol = 2),
    mu_beta = c(0,0),
    sigma_beta = c(1,1),
    mu_bg = c(-8.25, 0.066),
    sigma_bg = c(0.01, 0.01),
    a_cf = 1,
    b_cf = 1#,
    # h_bg = tx_dat[[tx_name]]$PFS_rate
  )

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

