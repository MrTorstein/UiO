from numpy import linspace, zeros_like, sin, sqrt
from matplotlib.pyplot import plot, title, xlabel, ylabel, legend, show

"""
Funksjonen er:
y(x, t) = A sin( k x - w(k) t )
der w(k) = c sqrt( k**2 + ( m c / h_ )**2 )
"""

c	= 1
h_	= 1
A	= 1
m	= 1
k_1	= 0.6
k_2	= 0.7
T	= 10
N	= 1e4
dt	= T / N
w	= lambda k: c * sqrt( k**2 + ( m * c / h_ )**2 )

def EC(w, k, T, N):
	t = linspace(0, T, N)
	x = zeros_like(t)
	y = zeros_like(t); y[0] = A * sin( k * x[0] - w(k) * t[0] )
	v = zeros_like(t)
	
	for i in xrange(len(t) - 1):
		a = - A * w(k)**2 * sin( k * x[i] - w(k) * t[i] )
		v[i + 1] = v[i] + a * dt
		x[i + 1] = x[i] + v[i + 1] * dt
		y[i + 1] = A * sin( k * x[i + 1] - w(k) * t[i] )
	
	return t, x, y, v

t1, x1, y1, v1 = EC(w, k_1, T, N)
t2, x2, y2, v2 = EC(w, k_2, T, N)

plot(x1, y1, x2, y2, x1, y1 + y2)
show()