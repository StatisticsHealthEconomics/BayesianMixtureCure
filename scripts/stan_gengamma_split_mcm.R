
# run stan mixture cure model
# generalised gamma distribution
# with separate censored and uncensored likelihood components


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
    N = nrow(tx_dat[[tx_name]]),
    n_unc = sum(tx_dat[[tx_name]]$pfs_event),
    n_cens = sum(1 - tx_dat[[tx_name]]$pfs_event),
    t_unc = tx_dat[[tx_name]]$pfs[tx_dat[[tx_name]]$pfs_event == 1],
    t_cens = tx_dat[[tx_name]]$pfs[tx_dat[[tx_name]]$pfs_event == 0],
    unc_idx = which(tx_dat[[tx_name]]$pfs_event == 1),
    cens_idx = which(tx_dat[[tx_name]]$pfs_event == 0),
    H = 2,
    X = matrix(c(rep(1, nrow(tx_dat[[tx_name]])),
               scale(tx_dat[[tx_name]]$PFSage)),
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
    h_bg = tx_dat[[tx_name]]$PFS_rate/12)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter <- 500#0
n_warmup <- 100#0

stan_base <-
  stan(
    file = here::here("stan", "gengamma_split_mixture_cure_model.stan"),
    data = data_list,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 12,
                   stepsize = 0.5),
    warmup = n_warmup,
    iter = n_iter,
    chains = 1)

stan_base

