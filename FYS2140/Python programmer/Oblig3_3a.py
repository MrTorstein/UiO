from numpy import linspace, zeros_like, pi, cos, sin, sqrt
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, legend, savefig, show

"""
Funksjonen er:
y(x, t = 1) = A sin( k x - w(k) )
der w(k) = c sqrt( k**2 + ( m c / h_ )**2 )
"""

c	= 1
h_	= 1
A	= 1
m	= 1
t	= 1
k_1	= 0.6
k_2	= 0.7
N	= 1e4
w	= lambda k: c * sqrt( k**2 + ( m * c / h_ )**2 )

x = linspace(- 100, 100, N)
y1 = A * sin(k_1 * x - w(k_1) * t)
y2 = A * sin(k_2 * x - w(k_2) * t)

figure()
plot(x, y1 + y2, x, y1, x, y2)
title("Superponering av to bolger")
xlabel("x")
ylabel("y")
legend(["Sammenlagt bolge", "y1(x, 1) = sin(0.6 * x - 0.6)", "y1(x, 1) = sin(0.7 * x - 0.7)"])
savefig("Figur_01.png")

figure()
plot(x, y1 + y2)
title("Superponering av to bolger")
xlabel("x")
ylabel("y")
legend(["y(x, 1) = sin(0.6 * x - 0.6) + sin(0.7 * x - 0.7)"])
savefig("Figur_02.png")

show()