#Exercise 20: 
# In this exercise we will study data from an experiment where one wanted to assess the toxicity of the substance rotenone. 
# Groups of about 50 insects were exposed to various doses of rotenone, and the number of insects that died at each dose level was recorded. 
# The data are available at the course web-page.

# The variables in the data set are coded as follows:
# - LOGDOSE     Logarithm of the concentration of rotenone (base 10 logarithm)
# - NUMBER      Number of insects in total
# - DEAD        Number of insects that died


# a) Compute the proportion of insects that died at each dose level and make a plot of the proportions versus dose of rotenone.

insects = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/insects.txt", header = T, na.strings = ".")
p = insects$DEAD / insects$NUMBER
plot(p ~ insects$LOGDOSE, xlab = "Log(dose)")
print("-----------------------------------------------------------")


# b) Fit a suitable regression model for the relation between the proportions of insects that died and the doses of rotenone. 

insects.fit = glm(cbind(DEAD, NUMBER - DEAD) ~ LOGDOSE, data = insects, family = binomial)
summary(insects.fit)

# Give the reasons for your choice of regression model and interpret the fitted model.

# I choose a logistic reg model since a model predicting values lower than 0 or higher than the amount insects in the group are impossible. 
# But also since the datapoints are quit linear.
print("-----------------------------------------------------------")


# c) Assess the fit of the model by including the fitted proportions on the plot from question a.
x = seq(0, 2, 0.01)
y = predict(insects.fit, list(LOGDOSE = x), type = "response")
lines(x,  y)

# Also give a formal goodness-of-fit test. 

anova(insects.fit, test = "Chisq")

# 
print("-----------------------------------------------------------")

# d) Use the fitted model to estimate LD50, i.e. the dose required to kill half the members of a tested population.

beta0 = insects.fit$coef[1]
beta1 = insects.fit$coef[2]
p = 0.5
LD50 = log(p / (1 - p)) - beta0 / beta1
print(paste("LD50 = ", LD50))
predict(insects.fit, list(LOGDOSE = LD50), type = "response")