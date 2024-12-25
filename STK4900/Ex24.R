library(survival)

melanoma = read.table("https://www.uio.no/studier/emner/matnat/math/STK4900/data/melanoma.txt", header = T)

melanoma.fit = survfit(Surv(lifetime, status == 1) ~ factor(sex), data = melanoma, conf.type = "plain")
plot(melanoma.fit, col = c("pink", "blue"), conf.int = T, ylim = c(0.43, 1.03))
legend(0.1, 0.5, c("Kvinne", "Mann"), lty = 1:1, col = c("pink", "blue"))

survdiff(Surv(lifetime, status == 1) ~ factor(sex), data = melanoma)

# With a p-value of a little over 0.01, it seems to be a small significant difference.
print("-------------------------------------------------------------------")


melanoma.fit = survfit(Surv(lifetime, status == 1) ~ factor(grthick), data = melanoma, conf.type = "plain")
plot(melanoma.fit, col = c("red", "green", "blue"), conf.int = T, ylim = c(0.2, 1.02))
legend(0.1, 0.31, c("0 - 1 mm", "1 - 2.5 mm", "5+ mm"), lty = 1:1, col = c("red", "green", "blue"))

survdiff(Surv(lifetime, status == 1) ~ factor(grthick), data = melanoma)

# This has a much larger significans.

print("-------------------------------------------------------------------")


melanoma.fit = survfit(Surv(lifetime, status == 1) ~ factor(ulcer), data = melanoma, conf.type = "plain")
plot(melanoma.fit, col = c("red", "green"), conf.int = T, ylim = c(0.2, 1.02))
legend(0.1, 0.31, c("present", "absent"), lty = 1:1, col = c("red", "green"))

survdiff(Surv(lifetime, status == 1) ~ factor(ulcer), data = melanoma)

# This does definantely have a much larger significans.
print("-------------------------------------------------------------------")


summary( coxph( Surv(lifetime, status == 1) ~ factor(sex), data = melanoma ) )
print("-------------------------------------------------------------------")

summary( coxph( Surv(lifetime, status == 1) ~ factor(ulcer), data = melanoma ) )
print("-------------------------------------------------------------------")

summary( coxph( Surv(lifetime, status == 1) ~ thickn, data = melanoma ) )
print("-------------------------------------------------------------------")

summary( coxph( Surv(lifetime, status == 1) ~ factor(grthick), data = melanoma ) )
print("-------------------------------------------------------------------")

summary( coxph( Surv(lifetime, status == 1) ~ logthick, data = melanoma ) )

# sex has small significans, just like earlier.
# ulcer has very large significans, just like earlier.
# all of the thick versions have quite large significans, but the numberical value thickn seems to have the smallest 95% confidence interval, so I would use that.

print("-------------------------------------------------------------------")

summary( coxph( Surv(lifetime, status == 1) ~ factor(sex) + factor(ulcer) + thickn + factor(grthick) + logthick, data = melanoma ) )