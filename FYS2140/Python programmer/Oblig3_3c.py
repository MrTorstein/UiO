from numpy import linspace, zeros_like, pi, cos, sin, sqrt, sum as Sum
from matplotlib.pyplot import plot, title, xlabel, ylabel, legend, savefig, show

"""
Funksjonen er:
y(x, t = 0) = A(k) cos( k x )
"""

c		= 1
h_		= 1
m		= 1
a 		= 1
t		= 0
k_abs	= 5
x_abs	= 20
n		= 1e4
N		= 1e4


k = linspace(- k_abs, k_abs, int(n))
A = a / (k**2 + a)
x = linspace(- x_abs, x_abs, int(N))
y = zeros_like(x)

for i in xrange(int(N)):
	y[i] = Sum(A * cos(k * x[i]))
	

plot(x, y)
title("Superponering av %d kosinus bolger" %n)
xlabel("x")
ylabel("y")
legend(["Sammenlagt bolge"])
#savefig("Figur_03.png")
show()