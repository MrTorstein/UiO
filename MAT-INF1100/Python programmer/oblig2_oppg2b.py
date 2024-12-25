from numpy import zeros, linspace
from matplotlib.pyplot import plot, show, figure, legend, savefig
from math import e

def f(t, x):
	return 1 - x**2

a, b = (0., 2.)
n = 5
dt = (b - a)/n

t = linspace(a, b, n + 1)
t2 = linspace(a, b, 1000)
x = zeros(n + 1)

for k in xrange(n):
	x[k + 1] = x[k] + dt*f(t[k], x[k])

xex = (e**(2*t2) - 1)/(e**(2*t2) + 1.)

plot(t, x, t2, xex)
legend(["Euler", "Eksakt"], loc="best")
savefig("oblig2oppg2b.png")
show()

"""
C:\Users\Torstein\Documents\UiO\MatInf1100\Python programmer>oblig2matinf_oppg2b.py

"""
