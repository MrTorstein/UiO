hers.sample = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v17/hers.sample.txt", header = T)

plot(hers.sample$age, hers.sample$sbp)
# The plot seem to tell us that there is almost no relation between age and blood pressure.

hers.fit.b = lm(sbp~age, data = hers.sample)
summary(hers.fit.b)
abline(hers.fit.b)

# The regression line rises to the right. It rises about 35/140*100 = 25% of its total value, which makes it likely that there is a corrolation between age and 
# blood pressure.
# The y-intercept is the blood pressure when a person is born, and the x-intercept is the age when the blood pressure is 0. Since this is a linear line and it rises
# 25% in about half a persons life, there are no positive age at which the blood pressure is zeros. This makes sense since no person has zero blood pressure.

print("____________________________________________________")

hers.fit.c = lm(sbp~I(age - 67), data = hers.sample)
summary(hers.fit.c)
abline(hers.fit.c)

# The least squares estimate is the same for both regressions.
# The interpretation of the intercept is the same as before.

print("____________________________________________________")

hers.fit.d = lm(sbp~I(age / 10), data = hers.sample)
summary(hers.fit.d)
abline(hers.fit.d)

# This least squares estimate is 10 times larger.
# I dont get a line from this for some reason.