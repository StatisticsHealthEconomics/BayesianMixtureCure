
# run stan mixture cure model
# generalised gamma distribution
# generating directly from priors


library(rstan)
library(shinystan)
library(dplyr)
library(ggplot2)


data_list <-
  list(
    a_cf = 0.8,
    b_cf = 7,
    a_mu = 0.3,
    b_mu = 0.6,
    a_Q = -0.2,
    b_Q = 0.6,
    a_scale = log(0.5),
    b_scale = 0.1,
    N_samples = 1)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

pred <-
  stan(
    file = here::here("stan", "gengamma_generated.stan"),
    data = data_list,
    iter = 1000, chains = 1,
    algorithm = "Fixed_param")

fit_stan <- extract(pred)


#######
# plot

fit_stan$S_pred <- array(fit_stan$S_pred, c(dim(fit_stan$S_pred), 1))

S_dat <- multimcm:::prep_S_data(fit_stan, tx_idx = 1)

ggplot(S_dat[[1]], aes(x = time, y = mean, group = type, colour = type)) +
  geom_line() +
  ylab("Survival") +
  ylim(0, 1) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = type),
              linetype = 0,
              alpha = 0.2) +
  theme_bw() +
  xlim(0,10)

