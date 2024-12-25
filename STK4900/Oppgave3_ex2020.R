path = "https://www.uio.no/studier/emner/matnat/math/STK4900/data/wheezing.txt"
wheezing = read.table(path, header = T)

#a)
n_subjects = length(wheezing$smoking)

n_smoke = sum(wheezing$smoking)
n_nosmoke = n_subjects - n_smoke

n_wheeze_nosmoke = sum(wheezing$wheeze[0:n_nosmoke])
n_wheeze_smoke = sum(wheezing$wheeze[(n_nosmoke + 1):n_subjects])

p_wheeze_nosmoke = n_wheeze_nosmoke / n_nosmoke
p_wheeze_smoke = n_wheeze_smoke / n_smoke

print(paste(p_wheeze_nosmoke, "of the children with non smoking mothers have wheezing"))
print(paste(p_wheeze_smoke, "of the children with smoking mothers have wheezing"))

p_wheeze_0 = (n_wheeze_nosmoke + n_wheeze_smoke) / n_subjects

se_wheeze_0 = sqrt((p_wheeze_0 - p_wheeze_0^2) / n_nosmoke + (p_wheeze_0 - p_wheeze_0^2) / n_smoke)

z = (p_wheeze_nosmoke - p_wheeze_smoke) / se_wheeze_0
print(paste("z =", abs(z)))

p_nowheeze_nosmoke = (n_nosmoke - n_wheeze_nosmoke) / n_nosmoke
p_nowheeze_smoke = (n_smoke - n_wheeze_smoke) / n_smoke

RR_nosmoke = p_wheeze_nosmoke / p_nowheeze_nosmoke
RR_smoke = p_wheeze_smoke / p_nowheeze_smoke

print(paste("Relative risk for non smoking =", RR_nosmoke))
print(paste("Relative risk for smoking =", RR_smoke))

OR_nosmoke = p_wheeze_nosmoke * (1 - p_nowheeze_nosmoke) / (p_nowheeze_nosmoke * (1 - p_wheeze_nosmoke))
OR_smoke = p_wheeze_smoke * (1 - p_nowheeze_smoke) / (p_nowheeze_smoke * (1 - p_wheeze_smoke))

print(paste("Odds ratio for non smoking =", OR_nosmoke))
print(paste("Odds ratio for smoking =", OR_smoke))


#b)
fit_logistic_smoke = glm(wheeze ~ smoking, data = wheezing, family = binomial)
summary(fit_logistic_smoke)

print(exp(unname(fit_logistic_smoke$coeff[2])))

fit_logistic_age = glm(wheeze ~ age, data = wheezing, family = binomial)
summary(fit_logistic_age)

fit_logistic_all = glm(wheeze ~ smoking + age, data = wheezing, family = binomial)
summary(fit_logistic_all)


#c)
library(gee)
path = "https://www.uio.no/studier/emner/matnat/math/STK4900/data/wheezingb.txt"
wheezing_b = read.table(path, header = T)

fit_logistic_all_b = gee(wheeze ~ smoking + age, id = id, data = wheezing_b, family = binomial, corstr = "unstructured")
summary(fit_logistic_all_b)