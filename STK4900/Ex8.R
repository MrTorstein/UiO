N_disp = c(0, 0, 1, 1, 2, 2, 4, 4, 5, 5, 6, 6, 7, 7)
sale = c(508.1, 498.4, 568.2, 577.3, 651.7, 657.0, 755.3, 758.9, 787.6, 792.1, 841.4, 831.8, 854.7, 871.4)

plot(N_disp, sale)

fit.line = lm(sale ~ N_disp)
summary(fit.line)
abline(fit.line, col = "blue")

fit.poly = lm(sale ~ poly(N_disp, 2))
predicted.intervals = predict(fit.poly, data.frame(N_disp = 0:7))

lines(0:7, predicted.intervals, col = "red")

# It is obvious that the second order line is a better fitt than the line. So this is the best model for this analesis.