library(MASS)

n = 400
rho = 0.90
m = matrix(c(0, 0), nrow = 2)
S = matrix(c(1, rho, rho, 1), nrow = 2)
obs = mvrnorm(n, m, S)
x = obs[, 1]
y = obs[, 2]
cor(x, y)
plot(x, y)

# The pearson corrolation is changing between 0.6 and -0.1. And is larger the closer the points are to a line from bottom left corner to top right corner.

# Same effect, but the values are mostly between 0.4 and 0.8.

# Interval is 0.8 to 1.0

# From this we see that the closer the rho is to 1.0, the smaller the 

# n = 100, rho = 0.30 => p_corr = (0.0, 0.5)

# n = 100, rho = 0.60 => p_corr = (0.5, 0.7)

# n = 100, rho = 0.90 => p_corr = (0.84, 0.92)

# n = 400, rho = 0.30 => p_corr = (0.25, 0.39)

# n = 400, rho = 0.60 => p_corr = (0.54, 0.68)

# n = 400, rho = 0.90 => p_corr = (0.89, 0.92)

# The larger the sample size, the smaler interval.