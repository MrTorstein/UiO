cancer_data = read.table("https://www.uio.no/studier/emner/matnat/math/STK4900/data/cancer.txt")
names(cancer_data) = c("age", "cig", "yl", "cancer")

cancer.fit_1 = glm(cancer ~ offset(log(yl)) + cig + age, data = cancer_data, family = poisson)
summary(cancer.fit_1)

# Very significant

expcoef = function(glmobj)
{
regtab  = summary(glmobj)$coef
expcoef = exp(regtab[, 1])

lower   = expcoef * exp(-1.96 * regtab[, 2])
upper   = expcoef * exp(1.96 * regtab[, 2])

cbind(expcoef, lower, upper)
}

expcoef(cancer.fit_1)

print("-----------------------------------------------------------------")

cancer.fit_2 = glm(cancer ~ offset(log(yl)) + cig + age + I(cig^2), data = cancer_data, family = poisson)
cancer.fit_3 = glm(cancer ~ offset(log(yl)) + cig + age + I(cig^2) + I(age^2), data = cancer_data, family = poisson)
cancer.fit_4 = glm(cancer ~ offset(log(yl)) + cig + age + I(cig^2) + I(age^2) + cig:age, data = cancer_data, family = poisson)
anova(cancer.fit_1, cancer.fit_2, cancer.fit_3, cancer.fit_4)

# The best model is model 3.