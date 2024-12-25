# Exercise 16:  confidence interval and tests for two proportions
# In this exercise we will use the results from the opinion polls from Norstat from February and March 2017 to investigate the change in the support for 
# some political parties.

# In February 2017 Norstat asked n1 = 935 individuals which party they would support if there had been election to the parliament tomorrow. 
# Of these y1 = 230 would have voted Høyre.
# One month later, in March, n2 = 927 persons were interviewed and y2 = 208 of these would have voted Høyre. See http://www.pollofpolls.no/?cmd=Maling&gallupid=3096
 

# Question a)
# We start out by estimating the change in the support for Høyre with a 95 % confidence interval  (cf. slide 6)
# Try to program such an interval yourself in R. A suggested solution is given at the bottom of this page. Comment on the results. 

y1  = 230
n1  = 935
p1  = y1 / n1
se1 = sqrt(p1 * (1 - p1) / n1)

y2  = 208
n2  = 927
p2  = y2 / n2
se2 = sqrt(p2 * (1 - p2) / n2)

p   = p1 - p2
se  = sqrt(se1^2 + se2^2)

mar = 1.96 * se
low = p - mar
upp = p + mar

cbind(p, mar, low, upp)

# The change is 2.16%, with a margin of error of +- 3.85%. This means the margin of error is bigger than the acual value of the change.


# Question b)
# We then test the null hypothesis that the support for Høyre has not changed from February to March (cf. slide 8)

p = (y1 + y2) / (n1 + n2)
se0 = sqrt(p * (1 - p) / n1 + p * (1 - p) / n2)
z = (p1 - p2) / se0
pval = 2 * (1 - pnorm(abs(z)))
cbind(z, pval)

# Perform these commands and comment on the results.
# Is the null hypothesis rejected or not? How does this relate to the confidence interval computed earlier?

# The p-value is 0.272, which is larger than the alpha value (= 5% = 0.05) so the null hypothesis can not be rejected.


# Question c)
# R has a command for comparing two proportions

hoyre = matrix(c(y1, y2, n1 - y1, n2 - y2), nrow = 2)   # give the data for Høyre in a 2x2 table (cf. slide 10)
prop.test(hoyre, correct = F)

# Perform these commands and check that the results agree with those obtained earlier.

# They agree.

# The prop.test-command give a chi squared statistic, not a z-value as we computed earlier. What is the relation between the two?
# chi^2 = z^2
 

# Question d)
# We will then take a look at the results for Senterpartiet (Sp). In February 80 of the 935  persons who were interviewed would have voted Senterpartiet; 
# while in March 101 of the 927 interviewed would have voted Senterpartiet.
# Estimating the change in the support for Senterpartiet with a 95 % confidence interval.

y3  = 80
n3  = 935
p3  = y3 / n3
se3 = sqrt(p3 * (1 - p3) / n3)

y4  = 101
n4  = 927
p4  = y4 / n4
se4 = sqrt(p4 * (1 - p4) / n4)

p   = p3 - p4
se  = sqrt(se3^2 + se4^2)

mar = 1.96 * se
low = p - mar
upp = p + mar

cbind(p, mar, low, upp)

# The change in persent is up 2.34%, and has the margin of error of 2.69%.

# Also test the null hypothesis that the support for Senterpartiet has not changed from February to March.

p = (y3 + y4) / (n3 + n4)
se0 = sqrt(p * (1 - p) / n3 + p * (1 - p) / n4)
z = (p3 - p4) / se0
pval = 2 * (1 - pnorm(abs(z)))
cbind(z, pval)

# What are your conclusions?

# The p-value is 0.088, which mean that we can not reject the null hypotesis, but the alternative is still significant.