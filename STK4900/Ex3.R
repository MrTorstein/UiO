time = c(28, -44, 29, 30, 24, 28, 37, 32, 36, 27, 26, 28, 29, 26, 27, 22, 23, 20, 25, 25, 36, 23, 31, 32, 24, 27, 33,
         16, 24, 29, 36, 21, 28, 26, 27, 27, 32, 25, 28, 24, 40, 21, 31, 32, 28, 26, 30, 27, 26, 24, 32, 29, 34, -2, 
         25, 19, 36, 29, 30, 22, 28, 33, 39, 25, 16, 23)

hist(time)
plot(ecdf(time))
boxplot(time)

print(paste("Mean value", mean(time)))
print(paste("Median", median(time)))

print(paste("Standard diviation", sd(time)))
quant = quantile(time)
Q1 = quant[2]
Q2 = quant[4]
IQR = Q2 - Q1
print(paste("Interquartile range", IQR))

conf = qt(0.95, 65)

print(paste("95% confidence interval", conf))
print(paste("Lower half", mean(time) - conf * (sd(time) / sqrt(66))))
print(paste("Upper half", mean(time) + conf * (sd(time) / sqrt(66))))

time_2 = time[time > (Q1 - 1.5 * IQR)]
time_2 = time_2[time_2 < (Q2 + 1.5 * IQR)]
n_2 = length(time_2)
conf_2 = qt(0.95, n_2 - 1)

print(paste("Lower half without outliers", mean(time_2) - conf_2 * (sd(time) / sqrt(n_2))))
print(paste("Upper half without outliers", mean(time[time_2]) + conf_2 * (sd(time) / sqrt(n_2))))

# The mean and median values of our dataset doesn't really match the true value very much, but the standard diviation is 
# quite large, and with that we do get the true value included. The confidence intervales do not contain the true value, 
# and are not very much better without the outliers.