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