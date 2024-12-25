#coding=utf-8
from numpy import linspace, zeros
from math import pi, e, sin
from matplotlib.pyplot import plot, show

t = linspace(-2*pi, 2*pi, 101)

x = zeros(101)
x[50] = 1 - e

dx = zeros(len(x))
for i in xrange(0, 50):
	dx[i + 50] = sin(i) - x[i + 50]*sin(i)
	x[i + 50 + 1] = dx[i + 50]*i+1 + x[i + 50]

#plot(t, x, t, dx)
#show()
#print len(x), len(t)
print dx, x