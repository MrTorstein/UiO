from numpy import random, sum as Sum, pi, sqrt, e, linspace
from matplotlib.pyplot import figure, plot, hist, title, legend, xlabel, ylabel, savefig, show

M = 10**4						# Antall mikrotilstander
N = 60							# Antall forskjellige spinn
d = 10**2						# Mellomrom mellom tilstander som blir plottet

liste = random.randint(-1, 1, (M, N)) * 2 + 1	# Lager en M X N liste av 1 eller -1
summer = Sum(liste, 1)							# Summerer spinnet til mikrotilstandene
spinn = linspace(-30, 30, M)
Omg = 2**N * sqrt(2 / (pi * N)) * e**(-spinn**2 / (2 * N))

figure()
plot(spinn, Omg * 10**(-14))
hist(summer, bins = "auto")
title("Histogramet fra tidligere med plott av multiplisiteten til ett spinn")
xlabel("Spinn")
ylabel("Multiplisitet / $10^{14}$ eller antall mikrotilstander")
legend(["$\Omega(S) \cdot 10^{-14}$", "Antall(S)"])
savefig("Figur104.png")

show()