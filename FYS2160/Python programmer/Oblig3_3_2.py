from numpy import linspace, zeros, exp, sqrt, pi
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, savefig, show

T1	= 300							# K
T2	= 600							# K
k	= 1.38 #* 10 ** (- 23)			# J / K
An	= 6.02 #* 10 ** (23)			# 1 / mol
m	= 0.028							# kg / mol

v	= linspace(0, 2000, 1E5)		# m/s
D1	= ( m / ( 2 * pi * k * T1 * An ) ) ** ( 3 / 2 ) * 4 * pi * v ** 2 * exp( - m * v ** 2 / ( 2 * k * T1 * An ) )
D2	= ( m / ( 2 * pi * k * T2 * An ) ) ** ( 3 / 2 ) * 4 * pi * v ** 2 * exp( - m * v ** 2 / ( 2 * k * T2 * An ) )

figure()
plot(v, D1 * 100, v, D2 * 100)
title("Maxwell distribusjons fuksjonen for nitrogen gass")
legend(["D(v) med T = %i K"%T1, "D(v) med T = %i K"%T2])
xlabel("Hastighet v [m/s]")
ylabel("D(v) [1 / v] * 100")
savefig("Figur03.png")
show()