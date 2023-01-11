
# kaplan-meier interval censoring

# test
tcut <- 3

time1 <- data_list$t
time1[time1 < tcut] <- 0

time2 <- data_list$t

event <- data_list$d
event[time1 < tcut & event == 1] <- 3

surv_dat <- Surv(time1, time2, event = event, type = "interval")

km <- survfit(surv_dat ~ 1)
plot(km)
