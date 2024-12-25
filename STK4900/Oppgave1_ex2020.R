path = "https://www.uio.no/studier/emner/matnat/math/STK4900/data/uffi.txt"
uffi = read.table(path, header = T)

#a)
CH2O_concentration = uffi$CH2O
UFFI_groups = uffi$UFFI

plot(UFFI_groups, CH2O_concentration)

t.test(uffi$CH2O[13:24], uffi$CH2O[0:12], var.equal = T)

line_fit_UFFI = lm(CH2O ~ UFFI, data = uffi)
summary(line_fit_UFFI)
abline(line_fit_UFFI)


#b)
Air_tightness = uffi$AIR

plot(Air_tightness, CH2O_concentration)

print(paste("Pearson correlation =", cor(Air_tightness, CH2O_concentration)))

line_fit_AIR = lm(CH2O ~ AIR, data = uffi)
summary(line_fit_AIR)
abline(line_fit_AIR)


#c)
line_fit_all = lm(CH2O ~ AIR + UFFI, data = uffi)
summary(line_fit_all)