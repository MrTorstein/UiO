from numpy import random, sum as Sum, pi, sqrt, e, linspace, log as ln
from matplotlib.pyplot import figure, plot, hist, title, legend, xlabel, ylabel, savefig, show

k = 8.62 * 10**(-5)						# Boltzmanns konstant [eV/K]
N = 60									# Antall forskjellige spinn

spinn = linspace(-N / 2, N / 2, 10**4)	# Array av spinn fra -30 til 30
SB = k * (N * ln(2) * 0.5 * ln(2 / (pi * N)) - spinn**2 / (2 * N))

figure()
plot(spinn, SB * 100)
title("Entropien til en paramagnet som en funksjon av netto spinn")
xlabel("Spinn")
ylabel("Entropi [eV/K]")
legend(["$S_B(S, N = 60) \cdot 100$"])
savefig("Figur105.png")

show()