# Exercise 15:  confidence interval for proportions
#At the lectures we looked an example concerning the opinion poll from February 2017 (cf. slide 2 from the lectures)
#We will in this exercise consider this example further.


# Question a)
# Of the n = 935 persons who were interviewed by Norstat, y = sum(y_i) = 309 would have voted Ap
# The following calculations reproduce the result from the lectures (cf. slide 4)

n = 935
y = 309
p_AP = y / n
se = sqrt(p_AP * (1 - p_AP) / n)
margin_AP = 1.96 * se
lower = p_AP - margin_AP
upper = p_AP + margin_AP
cbind(p_AP, margin_AP, lower, upper)

# Do the calculations and check that you get the result from the lectures.


# Question b)
# In the opinion poll, 122 of the persons interviewed would have voted Fremskrittspartiet (FrP) and 80 would have voted Senterpartiet (Sp).

# Repeat the calculations above for Fremskrittspartiet and Senterpartiet.
y_FRP = 122
p_FRP = y_FRP / n
se = sqrt(p_FRP * (1 - p_FRP) / n)
margin_FRP = 1.96 * se
lower = p_FRP - margin_FRP
upper = p_FRP + margin_FRP
cbind(p_FRP, margin_FRP, lower, upper)

y_SP = 80
p_SP = y_SP / n
se = sqrt(p_SP * (1 - p_SP) / n)
margin_SP = 1.96 * se
lower = p_SP - margin_SP
upper = p_SP + margin_SP
cbind(p_SP, margin_SP, lower, upper)

# How is the "margin of error" for these parties compared to the "margin of error" for Ap (cf. slide 4)?
print(margin_AP / p_AP * 100)
print(margin_FRP / p_FRP * 100)
print(margin_SP / p_SP * 100)

# The margin of error is much bigger for FRP and SP than it is for AP compared to the actual value of the estimated persentage. However, the margin it self is smaller.
# This means the margin of error sinks slower than the actual value does.