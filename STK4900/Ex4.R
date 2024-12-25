solvents = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/solvents.txt", header = T)

boxplot(rate~type, data = solvents)
# The boxplots tell me that the sorption rate for chemical one and two are higher than number three, and that they are probably close to eachother in value.

solvents$type = factor(solvents$type)
aov.solvents = aov(rate~type, data = solvents)
summary(aov.solvents)
# It seems that since the F value is quite large it is some difference in the average value for the sorption rates.