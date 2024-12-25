pef = c(494, 395, 516, 434, 476, 413, 442, 433)
mpef = c(512, 430, 520, 428, 500, 364, 380, 445)
person = 1:8

plot(person, pef)
plot(person, mpef)

PC = cor(pef, mpef)
print(PC)

# Pearson corrolation = 0.8154886

cor.test(pef, mpef)

plot(pef, mpef, pch = 19)

fit = lm(mpef ~ pef)
summary(fit)
abline(fit)

# Least sqares fit of pef is 1.1642.

print(1.1642 * sd(pef) / mean(pef))

# Res = 0.106893