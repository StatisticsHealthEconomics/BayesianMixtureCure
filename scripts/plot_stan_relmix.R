
# stan output plots
# relative survival mixture model


library(ggplot2)
library(dplyr)
library(reshape2)

fit_stan <- extract(stan_base)

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

