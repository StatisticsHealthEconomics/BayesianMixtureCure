
# run stan mixture cure model
# generalised gamma distribution


library(rstan)
library(shinystan)
library(dplyr)

load("../bgfscure/data/surv_input_data.RData")

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
    mu_beta = c(0.3, 0),
    sigma_beta = c(0.6, 0.1),
    a_Q = -0.2,
    b_Q = 0.6,
    a_scale = log(1.5),
    b_scale = 0.1,
    a_cf = 3,
    b_cf = 12,
    ## background
    # mu_bg = c(-8.25, 0.066),
    # sigma_bg = c(0.01, 0.01),
    h_bg = tx_dat[[tx_name]]$PFS_rate/12
  )

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter <- 500#0
n_warmup <- 100#0

stan_base <-
  stan(
    file = here::here("stan", "gengamma_mixture_cure_model.stan"),
    data = data_list,
    warmup = n_warmup,
    iter = n_iter,
    chains = 1)

stan_base

