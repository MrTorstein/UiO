from numpy import linspace, sqrt, pi, e, sinh
from matplotlib.pyplot import figure, subplot, plot, legend, title, xlabel, ylabel, axis, savefig, show

V_0	= 34		# [MeV]
dx	= 17		# [fm]
m	= 3727		# 3.727 [GeV / c^2] = 3.727 * 10^3 [MeV / c^2] = 3727 [MeV / c^2]
h_c	= 1.973e2	# 0.1973 [eV um] = 0.1973 * 10^9 * 10^-6 [MeV fm] = 1.973 * 10^2 [MeV fm]

E	= linspace(1e-10, 34, 10000)
T	= 1. / ( 1 + ( V_0 ** 2 ) / ( 4 * E * ( V_0 - E ) ) * sinh( dx * sqrt( 2 * m * ( V_0 - E ) ) / ( h_c ) ) ** 2 )

figure()
s1 = subplot(211)
plot(E, T)

legend(["T(E)"])
title("Plott av tran.koeffisient mot energi til alphapartikkelen")
ylabel("T [Enhetsloos]")

s2 = subplot(212)
plot(E, T)
legend(["T(E)"])
xlabel("E [MeV]")
ylabel("T [Enhetsloos]")
axis((33.5, 34, - 0.0001, 0.0021))
savefig("Figur_07.png")


show()