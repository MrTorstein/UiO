crabs = read.table("https://www.uio.no/studier/emner/matnat/math/STK4900/data/crabs.txt", header = TRUE)
summary(crabs)

fit.crabs_width = glm(y ~ width, data = crabs, family = binomial)
summary(fit.crabs_width)
print("-------------------------------------------------------")

OR = exp(fit.crabs_width$coef[2] * 1)

print("-------------------------------------------------------")

expcoef = function(glmobj)
{
regtab = summary(glmobj)$coef
expcoef = exp(regtab[, 1])
lower = expcoef * exp(-1.96 * regtab[, 2])
upper = expcoef * exp(1.96 * regtab[, 2])
cbind(expcoef, lower, upper)
}

expcoef(fit.crabs_width)

print("-------------------------------------------------------")

summary(glm(y ~ weight, data = crabs, family = binomial))
summary(glm(y ~ factor(color), data = crabs, family = binomial))
summary(glm(y ~ factor(spine), data = crabs, family = binomial))

print("-------------------------------------------------------")

fit.crabs_multi = glm(y ~ width + weight, data = crabs, family = binomial)
summary(fit.crabs_multi)
fit.crabs_corre = glm(y ~ width:weight, data = crabs, family = binomial)
summary(fit.crabs_corre)

print("-------------------------------------------------------")

fit.crabs_correlation = glm(y ~ weight + width + factor(spine) + factor(color) + weight:width + weight:factor(spine) + weight:factor(color) + width:factor(spine)
+ width:factor(color) + factor(spine):factor(color) + weight:width:factor(spine) + weight:width:factor(color) + weight:width:factor(spine):factor(color), 
data = crabs, family = binomial)
anova(fit.crabs_correlation, test = "Chisq")