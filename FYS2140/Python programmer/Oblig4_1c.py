from numpy import linspace, zeros_like, sqrt, exp, abs
from matplotlib.pyplot import plot, xlabel, ylabel, legend, title, savefig, show

"""
Variabler:
- lamb	= konstant
- A		= konstant = sqrt(lamb)
- t		= tidspunkt
- x		= posisjon i en dimensjon
- sig	= standaravvik
- N		= antall tidssteg
- psi	= |psi(x)|^2, der psi(x) = Ae^(- lamb * |x|) saa lenge t = 0
"""

lamb	= 2.
A		= sqrt(lamb)
t		= 0
abs_x	= 2.5
N		= 1e4
x		= linspace(- abs_x, abs_x, N)
sig		= 1. / ( sqrt(2) * lamb )
psi		= A**2 * exp(- 2 * lamb * abs(x))
forv_x	= 0

plot(x, psi, [forv_x - sig, forv_x - sig + 1e-10], [- 1e-1, 2], [forv_x + sig, forv_x + sig + 1e-10], [- 1e-1, 2])
xlabel("|psi|^2 [relv. enhet]")
ylabel("x [relv. enhet]")
legend(["A**2 e**(-2*lamb*abs(x))", "<x> - sigma", "<x> + sigma"])
title("Grafen til $|\psi(x)|^{2}$")
savefig("Figur_04.png")
show()