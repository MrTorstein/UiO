#Exercise 18:  a second look at logistic regression
# In this exercise we will use data from the Western Collaborative Group Study (WCGS), 
# a large epidemiological study designed to study risk factors for coronary heart disease (CHD). 
# More than 3000 men were recruited to the study, and a number of (potential) risk factors were recorded at entry. 
# The men were then followed for about ten years and it was recorded if they developed CHD or not. 
# The WCGS data are used in many of the examples in Section 3.4 and Chapter 6 of the course text book by Vittinghoff et al.

# You may read the WCGS data into R by the command:

wcgs = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/wcgs.txt", sep = "\t", header = T, na.strings = ".")

# In this exercise we will restrict ourselves to study the effect of smoking and age on the risk for CHD.
 

# Question a)
# We first fit a logistic regression model with smoke as covariate:

fit.smoke = glm(chd69~smoke, data = wcgs, family = binomial)
summary(fit.smoke)

# Is there a significant effect of smoking on the risk of developing CHD? (Cf. slide 11)

# Yes, the p-value is very low, which indicates a significans

# Compute a 95% confidence interval for the regression coefficient for smoking (cf slide 9)

beta_hat    = fit.smoke$coefficient[2]
se          = summary(fit.smoke)$coef[2, 2]
margin      = 1.96 * se
lower       = beta_hat - margin
upper       = beta_hat + margin
cbind(beta_hat, margin, lower, upper)

# Also estimate the odds ratio for smoking and determine a 95% confidence interval for the odds ratio (cf slide 9).

OR          = exp(beta_hat)
lower_OR    = exp(lower)
upper_OR    = exp(upper)
cbind(OR, lower_OR, upper_OR)
print("----------------------------------------------------------------")

# Question b)
# In order to estimate the odds ratio with 95% confidence limits from a fitted logistic model, you may use the function (cf slide 10)

expcoef = function(glmobj)
{
regtab  = summary(glmobj)$coef
expcoef = exp(regtab[, 1])

lower   = expcoef * exp(-1.96 * regtab[, 2])
upper   = expcoef * exp(1.96 * regtab[, 2])

cbind(expcoef, lower, upper)
}

# Use this function to estimate the odds ratio for smoking with 95% confidence limits.

expcoef(fit.smoke)

# Check that you get the same results as in question a.

# I did!
print("----------------------------------------------------------------")


# Question c)
# We then use logistic regression to study the effect of age (at entry to the study) for the risk of developing CHD:

fit.age = glm(chd69~age, data = wcgs, family = binomial)
summary(fit.age)

# Use the expcoef-function to estimate the odds ratio for a one-year increase in age with 95% confidence interval.

expcoef(fit.age)
print("----------------------------------------------------------------")

# Question d)
# When interpreting the effect of age, it may be reasonable to give the odds ratio corresponding to a ten-year increase in age (rather than a one-year increase as we did in question c). The easiest way to achieve this is to fit the logistic model using age/10 as covariate:

fit.age10 = glm(chd69~I(age / 10), data = wcgs, family = binomial)
summary(fit.age10)

# Compare the estimates and standard errors of this model with the one in question c. What do you see?

# The estimated value and standard error is just 10 times larger. z-value and p-value is not affected.

# You may then use the expcoef-function to estimate the odds ratio for a ten-year increase in age with 95% confidence interval. (Why?)

expcoef(fit.age10)

# The expcoef funk takes a glm object, so you have to recalculate the fitted model first?