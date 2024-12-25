library(survival)

cirrhosis = read.table("https://www.uio.no/studier/emner/matnat/math/STK4900/data/cirrhosis.txt", header = T)

cirrhosis.fit_treat = survfit(Surv(time, status) ~ treat, data = cirrhosis, conf.type = "plain")
plot(cirrhosis.fit_treat, col = c("green", "dark green"), conf.int = T)
legend(3500, 1, c("Prednisone", "Placebo"), lty = 1:1, col = c("green", "dark green"))

cirrhosis.fit_sex = survfit(Surv(time, status) ~ sex, data = cirrhosis, conf.type = "plain")
plot(cirrhosis.fit_sex, col = c("pink", "light blue"), conf.int = T)
legend(3500, 1, c("Female", "Male"), lty = 1:1, col = c("pink", "light blue"))

cirrhosis.fit_ascite = survfit(Surv(time, status) ~ asc, data = cirrhosis, conf.type = "plain")
plot(cirrhosis.fit_ascite, col = c("brown", "yellow", "orange"), conf.int = T)
legend(3000, 1, c("No lvl. of ascite", "Slight lvl. of ascite", "Marked lvl. of ascite"), lty = 1:1, col = c("brown", "yellow", "orange"))

cirrhosis.fit_age = survfit(Surv(time, status) ~ agegr, data = cirrhosis, conf.type = "plain")
plot(cirrhosis.fit_age, col = c("red", "purple", "blue"), conf.int = T)
legend(3000, 1, c("< 50 years old", "50 - 65 years old", "> 65 years old"), lty = 1:1, col = c("red", "purple", "blue"))

print("-------------------------------------------------------")

survdiff(Surv(time, status) ~ treat,    data = cirrhosis)
survdiff(Surv(time, status) ~ sex,      data = cirrhosis)
survdiff(Surv(time, status) ~ asc,      data = cirrhosis)
survdiff(Surv(time, status) ~ agegr,    data = cirrhosis)

print("-------------------------------------------------------")

summary(coxph(Surv(time, status) ~ treat + sex + asc + age, data = cirrhosis))