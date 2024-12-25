# Read values from course webpage #
Blood = read.table(
"https://www.uio.no/studier/emner/matnat/math/STK4900/v20/mandatory/blood.txt", header = T)

# Make two boxplots #
boxplot(Blood$Bloodpr)
boxplot(Blood$age)

# Print the values of the different quantiles #
quantile(Blood$Bloodpr)
quantile(Blood$age)

# Standard diviations #
print(paste("Standard diviation for Bloodpr", sd(Blood$Bloodpr)))
print(paste("Standard diviation for age", sd(Blood$age)))

