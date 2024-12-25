#At the lectures we looked an example on bone mineral density (cf. slide 25 from the lectures)
#We will in this exercise see how the computations for this example may be done in R.


#Start by reading the data into R. This may be done by the command:
cont = c(0.228, 0.207, 0.234, 0.220, 0.217, 0.228, 0.209, 0.221, 0.204, 0.220, 0.203, 0.219, 0.218, 0.245, 0.210)
test = c(0.250, 0.237, 0.217, 0.206, 0.247, 0.228, 0.245, 0.232, 0.267, 0.261, 0.221, 0.219, 0.232, 0.209, 0.255)


# Find the means and standard deviations, and check that you get the same results as in the lectures (slide 25)
print(paste("Mean control", mean(cont)))
print(paste("Standard diviation control", sd(cont)))
print(paste("Mean test", mean(test)))
print(paste("Standard diviation test", sd(test)))

# Use the command "t.test" to compute the confidence interval, t-statistic and P-value:

t.test(test, cont , var.equal = T)

# Make sure that you understand the output from the command and check that you get the same results as in the lectures (slides 27 and 28)

# Optional: Use the formulas given on slide 26 to compute the pooled standard deviation, the standard error of the effect of treatment, and the 95% confidence interval.

# Optional: Use the formulas given on slide 28 to compute the t-statistic and the P-value. 