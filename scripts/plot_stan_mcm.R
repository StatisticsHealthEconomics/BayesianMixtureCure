
# stan output plots
# relative survival mixture model


library(ggplot2)
library(dplyr)
library(reshape2)

fit_stan <- extract(stan_base)


# plots

plot(hist(fit_stan$curefrac,  breaks = 50))
abline(v = mean(fit_stan$curefrac), col = "red")

# survival curves

mean_rate0 <- mean(fit_stan$rate0)

time_horizon <- 6

dat <-
  data.frame(time = seq(0, time_horizon, length.out=100)) %>%
  mutate(rate0 = mean_rate0*time,
         px0 = 1 - pexp(time, rate = rate0)) %>%
  melt(id.vars = "time",
       measure.vars = "px0")

p1 <-
  ggplot(dat, aes(x=time, y=value, group=variable)) +
  geom_line()

p1

# credible intervals
q_rate0 <- quantile(fit_stan$rate0, probs = c(0.025, 0.975))

ribbon_df <-
  data.frame(time = seq(0, time_horizon, length.out=100)) %>%
  mutate(lower0 = q_rate0[1]*time,
         upper0 = q_rate0[2]*time,
         lS0 = 1 - pexp(time, rate = lower0),
         uS0 = 1 - pexp(time, rate = upper0))

p1 <-
  p1 +
  geom_ribbon(data = ribbon_df, aes(x = time, ymin = lS0, ymax = uS0),
              inherit.aes = FALSE,
              fill = "lightgrey",
              alpha = 0.5)

p1


#################
# using multimcm

library(ggplot2)
library(multimcm)

fit_stan <- extract(stan_base)

# extend dimension for a single treatment
fit_stan$S_pred <- array(fit_stan$S_pred, c(dim(fit_stan$S_pred), 1))

S_dat <- multimcm:::prep_S_data(fit_stan, tx_idx = 1)

gg <-
  ggplot(S_dat[[1]], aes(x = time, y = mean, group = type, colour = type)) +
  geom_line() +
  ylab("Survival") +
  ylim(0, 1) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = type),
              linetype = 0,
              alpha = 0.2) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_bw()

## kaplan-meier
library(survival)

km <- survfit(Surv(data_list$t, data_list$d) ~ 1)
# km <- survfit(Surv(c(data_list$t_unc, data_list$t_cens), c(rep(1,data_list$n_unc), rep(0,data_list$n_cens))) ~ 1)
# km <- survfit(Surv(dat$t, dat$d) ~ 1) |> plot()
km_data <- data.frame(surv = km$surv,
                      time = km$time)

gg + geom_step(aes(x = time, y = surv),
          linewidth = 1,
          data = km_data,
          inherit.aes = FALSE) +
  xlim(0,70)


