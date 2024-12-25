from numpy import zeros, linspace
from matplotlib.pyplot import plot, show, figure, legend, savefig
from math import e

def f(t, x):
	return 1 - x**2

a, b = (0., 2.)
n = 5
dt = (b - a)/n

t = linspace(a, b, n + 1)
tex = linspace(a, b, 1000)
x = zeros(n + 1)

for k in xrange(n):
	x[k + 1] = x[k] + dt*(1 - x[k]**2)

xex = (e**(2*tex) - 1)/(e**(2*tex) + 1.)

xm = zeros(n + 1)

for k in xrange(n):
	x_1 = xm[k] + dt*f(t[k], xm[k])/2.
	xm[k + 1] = xm[k] + dt*f(t[k], x_1)

plot(t, x, tex, xex, t, xm)
legend(["Euler", "Eksakt", "EulerMid"], loc="best")
savefig("oblig2oppg2c.png")
show()

"""
C:\Users\Torstein\Documents\UiO\MatInf1100\Python programmer>oblig2matinf_oppg2c.py

"""