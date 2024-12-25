olympic = read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/data/olympic.txt", sep = "\t", header = TRUE)

fit.olympic = glm(Total2000 ~ offset(Log.athletes) + Total1996 + Log.population + GDP.per.cap, data = olympic, family = poisson)
summary(fit.olympic)

print("-------------------------------------------------------")

fit._1996 = glm(Total2000 ~ offset(Log.athletes) + Log.population + GDP.per.cap, data = olympic, family = poisson)
summary(fit._1996)