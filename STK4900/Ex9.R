Income = c(47.35, 29.26, 52.14, 32.15, 40.86, 19.19, 27.23, 25.60, 54.14, 26.72, 38.84, 32.99, 32.95, 21.69, 27.90, 56.70, 37.69, 39.94)    # x-value
Risk = c(7, 5, 10, 6, 4, 5, 4, 6, 9, 5, 2, 7, 4, 3, 5, 1, 8, 6) # y-value
Insurance = c(140, 45, 180, 160, 90, 10, 35, 35, 190, 35, 75, 70, 55, 10, 40, 175, 95, 95)  # y-value
xval = seq(15, 60, by = (60 - 15) / 17)

plot(Income, Risk)  # Scatterplotting

Line.fit.risk = lm(Risk ~ Income)   # Linear regression
abline(Line.fit.risk, col = "blue")

Poly1.fit.risk = lm(Risk ~ poly(Income, 2))  # 2nd degree polyfit
Predicted = predict(Poly1.fit.risk, data.frame(Income = xval))
lines(xval, Predicted, col = "red")

Poly2.fit.risk = lm(Risk ~ poly(Income, 3))  # 3rd degree polyfit
predicted = predict(Poly2.fit.risk, data.frame(Income = xval))
lines(xval, predicted, col = "green")


plot(Income, Insurance) # Scatterplott

Line.fit.ins = lm(Insurance ~ Income)   # Linear regression
abline(Line.fit.ins, col = "blue")

Poly.fit.ins = lm(Insurance ~ poly(Income, 2))
predicted = predict(Poly.fit.ins, data.frame(Income = xval))
lines(xval, predicted, col = "red")

# It is at least obvious that the income and insurance are linearly dependant on eachother.
# As for The risk, it looks like it might be related by a third degree polynomial, but I would say it is not very likely.