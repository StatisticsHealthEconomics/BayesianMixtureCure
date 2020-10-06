
#
# stan output plots
#


library(ggplot2)
library(dplyr)
library(reshape2)

fit_stan <- extract(stan_base)

plot(hist(fit_stan$curefrac,  breaks = 50))
abline(v = mean(fit_stan$curefrac), col = "red")

# survival curves

mean_rate0 <- mean(fit_stan$rate0)
mean_rate1 <- mean(fit_stan$rate1)

time_horizon <- 2
# time_horizon <- 60

dat <-
  data.frame(time = seq(0, time_horizon, length.out=100)) %>%
  mutate(rate0 = mean_rate0*time,
         rate1 = mean_rate1*time,
         px0 = 1 - pexp(time, rate = rate0),
         px1 = 1 - pexp(time, rate = rate1)) %>%
  melt(id.vars = "time",
       measure.vars = c("px0", "px1"))

p1 <-
  ggplot(dat, aes(x=time, y=value, group=variable)) +
  geom_line()

p1

# credible intervals
q_rate0 <- quantile(fit_stan$rate0, probs = c(0.025, 0.975))
q_rate1 <- quantile(fit_stan$rate1, probs = c(0.025, 0.975))

ribbon_df <-
  data.frame(time = seq(0, 2, length.out=100)) %>%
  mutate(lower0 = q_rate0[1]*time,
         upper0 = q_rate0[2]*time,
         lower1 = q_rate1[1]*time,
         upper1 = q_rate1[2]*time,
         lS0 = 1 - pexp(time, rate = lower0),
         uS0 = 1 - pexp(time, rate = upper0),
         lS1 = 1 - pexp(time, rate = lower1),
         uS1 = 1 - pexp(time, rate = upper1))

p1 <-
  p1 +
  geom_ribbon(data = ribbon_df, aes(x = time, ymin = lS0, ymax = uS0),
              inherit.aes = FALSE,
              fill = "lightgrey",
              alpha = 0.5) +
  geom_ribbon(data = ribbon_df, aes(x = time, ymin = lS1, ymax = uS1),
              inherit.aes = FALSE,
              fill = "lightgrey",
              alpha = 0.5)

p1

## true values
# mean_rate0 <- 0.5
# mean_rate1 <- 2
# dat_true <-
#   data.frame(time = seq(0, 2, length.out=100)) %>%
#   mutate(rate0 = mean_rate0*time,
#          rate1 = mean_rate1*time,
#          px0 = 1 - pexp(time, rate = rate0),
#          px1 = 1 - pexp(time, rate = rate1)) %>%
#   melt(id.vars = "time",
#        measure.vars = c("px0", "px1"))
# 
# p1 +
#   geom_line(data = dat_true, aes(x=time, y=value, group=variable),
#             inherit.aes = FALSE, col = "red")

