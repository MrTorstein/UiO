# Exercise21: Fitting Poisson models

# The data in the table has all been taken from^10. The yeast cell relate to the number of blood cells counted by a so-called hemacytometer in small squares of blood.
# The moss shoots are similarly the number of a certain type of organisms in quadrats laid out by ecologists in some environments of interest.
# The azuki bean wevils are the number of wevils which have entered the beans, fed and purpated inside them, finally emerging and leaving behind holes. 
# The number of holes has been counted on the various beans of the experiment.
# In each of the three cases, fit the Poisson model, compute the "expected" tables and the coefficient of dispersion (CD).
# Judge whether the data seems to be Poisson distributed and speculate on reasons for deviations.

freq    = seq(0, 6, 1)
cell    = c(75, 103, 121, 54, 30, 13, 4)
moss    = c(100, 9, 6, 8, 1, 0, 2)
bean    = c(61, 50, 1)

n_cell = sum(cell)
y_cell = sum(cell * freq) / n_cell
s2_cell = sum(cell * (freq - y_cell)^2) / (n_cell - 1)
CD = s2_cell / y_cell
df = length(freq) - 2

print(paste("CD = ", CD))
print("________________________________________________")
p_cell_ = dpois(-1 : df + 1, y_cell)
print(p_cell_)
print("________________________________________________")
p_cell = c(dpois(0 : df, y_cell), 1 - ppois(df, y_cell))
print(p_cell)
print("________________________________________________")
E_cell = n_cell * p_cell_
cbind(freq, cell, E_cell)