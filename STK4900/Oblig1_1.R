# Read values from course webpage #
Values = read.table(
"https://www.uio.no/studier/emner/matnat/math/STK4900/v20/mandatory/no2.txt", header = T)

# Make two boxplots #
boxplot(Values$log.cars)
boxplot(Values$log.no2)

# Print the values of the different quantiles #
quantile(Values$log.cars)
quantile(Values$log.no2)

# Standard diviations #
print(paste("Standard diviation for log.cars", sd(Values$log.cars)))
print(paste("Standard diviation for log.no2", sd(Values$log.no2)))

# Scatterplot #
plot(log.no2~log.cars, data = Values)

# Making the linear regression model #
lrm = lm(log.no2~log.cars, data = Values)

# Printing a summary #
summary(lrm)

# Plotting the line in a new scatterplot #
plot(log.no2~log.cars, data = Values, col = "blue")
abline(lrm, col = "red")

# Producing and plotting the residual values, and a horizon #
res = resid(lrm)
plot(res)
abline(0, 0)

# Multiple linear regression #
mlrm = lm(log.no2~log.cars + temp + log(wind.speed) + hour.of.day, data = Values)
summary(mlrm)


# Trying to linearly plot something I believe is a 3D plot, and this gives nothing of value #
beta0 = mlrm$coefficients[1]
beta1 = mlrm$coefficients[2]
beta3 = mlrm$coefficients[4]

plot(log.no2~log.cars, data = Values, col = "blue")
abline(a = beta0, b = (beta1 + beta3), col = "red")

