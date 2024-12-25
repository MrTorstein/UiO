# You may copy the commands below from the web-browser into the command-window of R (or into a R-script)
# A line that starts with # is a comment, and R will disregard such lines.
# At the lectures we looked at an example with the age of mineral samples (cf. slide 3 from the lectures)
# We will in this exercise see how the computations for this examples may be done in R.


#Start by reading the data into R. This may be done by the command:

age = c(249, 254, 243, 268, 253, 269, 287, 241, 273, 306, 303, 280, 260, 256, 278, 344, 304, 283, 310)

 
#Compute mean, median and standard deviation:

mean(age)

median(age)

sd(age)

# Check that you get the same result as in the lectures (cf slide 3)


# Make a histogram (cf. slide 4)

hist(age)
 

# Plot the empirical distribution function (cf. slide 4)

plot(ecdf(age))                                         # Basic plot

plot(ecdf(age),verticals = T, do.points = F)          # Nicer looking plot


# Compute min, first quartile, median, third quartile, and max (cf. slide 5)

quantile(age)


# Make a boxplot (cf. slide 5)

boxplot(age)

 
# We will then consider confidence intervals and hypothesis testing.
# We will illustrate direct calculations of the quantities involved as well as the use of a special R-command.


# Compute the 97.5% percentile of the t-distribution with 18 degrees of freedom:

qt(0.975, 18)

# Compute lower and upper limit of the 95% confidence interval:

mean(age) - qt(0.975,18)*(sd(age)/sqrt(19))     # lower limit

mean(age) + qt(0.975,18)*(sd(age)/sqrt(19))     # upper limit


# Check that you get the same result as in the lectures (cf slide 18)

# Compute t-statistic:

tstat = (mean(age) - 265) / (sd(age) / sqrt(19))       #t-statistic

tstat

 

# Compute P-value:

1 - pt(tstat, 18)

# Check that you get the same result as in the lectures (cf slide 22)

# R has readymade commands for t-tests with corresponding confidence intervals.
# Use the command "t.test" to compute the confidence interval (this gives a two-sided test):

t.test(age, mu = 265)
 
# Use the command "t.test" to compute a one-sided test (this gives a one-sided confidence interval):

t.test(age, alternative = "greater", mu = 265)

# Check that you get the same results as above. 