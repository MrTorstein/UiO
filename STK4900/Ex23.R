library(survival)

gehan = read.table("https://www.uio.no/studier/emner/matnat/math/STK4900/data/gehan.txt", header = T, na.strings = ".")

group = c(rep(1, 21), rep(2, 21))
gehan_both.fit = survfit(Surv(time, cens) ~ treat, data = gehan, conf.type = "none")
summary(gehan_both.fit)
print(gehan_both.fit)


gehan_control.fit = survfit(Surv(time[1:21], cens[1:21]) ~ 1, data = gehan, conf.type = "plain")
plot(gehan_control.fit, col = "blue", lty = 1:1, xlim = c(0, 36))

gehan_drug.fit = survfit(Surv(time[22:length(gehan$time)], cens[22:length(gehan$time)]) ~ 1, data = gehan, conf.type = "plain")
lines(gehan_drug.fit, col = "red")

legend(1, 0.2, c("control", "drug"), col = c("blue", "red"), lty = 1:1, )

survdiff(Surv(time, cens) ~ group, data = gehan)

# a p-value of 4.17 * 10^-5 tells me it is a significant difference between the two groups.


summary( coxph( Surv(time, cens) ~ factor(treat), data = gehan ) )

# Seems that there is a significant use to the drug.