path = "https://www.uio.no/studier/emner/matnat/math/STK4900/data/drvisits.txt"
drv = read.table(path, header = T)

#a)
drv_fit = glm(numvisit ~ reform, data = drv, family = poisson)
summary(drv_fit)

coef  = summary(drv_fit)$coef
expcoef = exp(coef[, 1])
lower   = expcoef * exp(-1.96 * coef[, 2])
upper   = expcoef * exp(1.96 * coef[, 2])

cbind(expcoef, lower, upper)


#b)
drv_fit_all = glm(numvisit ~ age + educ + married + badh + loginc + reform, data = drv, family = poisson)
summary(drv_fit_all)


#c)
drv_fit_nonlin_0 = glm(numvisit ~ age + educ + married + badh + loginc + reform, data = drv, family = poisson)
drv_fit_nonlin_1 = glm(numvisit ~ age + educ + married + badh + loginc + reform + I(age^2), data = drv, family = poisson)
drv_fit_nonlin_2 = glm(numvisit ~ age + educ + married + badh + loginc + reform + I(age^2) + I(educ^2), data = drv, family = poisson)
drv_fit_nonlin_3 = glm(numvisit ~ age + educ + married + badh + loginc + reform + I(age^2) + I(educ^2) + I(loginc^2), data = drv, family = poisson)
drv_fit_nonlin_4 = glm(numvisit ~ age + educ + married + badh + loginc + reform + I(age^2) + I(educ^2) + I(loginc^2) + I(log(age)), data = drv, family = poisson)
drv_fit_nonlin_5 = glm(numvisit ~ age + educ + married + badh + loginc + reform + I(age^2) + I(educ^2) + I(loginc^2) + I(log(age)) + I(log(educ)), data = drv, family = poisson)
drv_fit_nonlin_6 = glm(numvisit ~ age + educ + married + badh + loginc + reform + I(age^2) + I(educ^2) + I(loginc^2) + I(log(age)) + I(log(educ)) + I(log(loginc)), data = drv, family = poisson)
anova(drv_fit_nonlin_0, drv_fit_nonlin_1, drv_fit_nonlin_2, drv_fit_nonlin_3, drv_fit_nonlin_4, drv_fit_nonlin_5, drv_fit_nonlin_6, test = "Chisq")


#d)
average = mean(drv$numvisit)

s2 = sum((drv$numvisit - average)^2) / (length(drv$numvisit) - 1)

print(paste("CD =", s2 / average))

drv_fit_quasi = glm(numvisit ~ age + educ + married + badh + loginc + reform, data = drv, family = quasipoisson)
summary(drv_fit_quasi)