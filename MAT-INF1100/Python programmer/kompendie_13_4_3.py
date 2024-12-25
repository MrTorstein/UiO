#coding=utf-8

#a)
from numpy import zeros, linspace, pi, e, sin
h1 = 1
h2 = 2
h3 = 5
h4 = 10
t1 = linspace(0, 2*pi, h1 + 1)
t2 = linspace(0, 2*pi, h2 + 1)
t3 = linspace(0, 2*pi, h3 + 1)
t4 = linspace(0, 2*pi, h4 + 1)

x1 = zeros(h1 + 1)
x1[0] = 2 + e
x2 = zeros(h2 + 1)
x2[0] = 2 + e
x3 = zeros(h3 + 1)
x3[0] = 2 + e
x4 = zeros(h4 + 1)
x4[0] = 2 + e

def f(t, x):
	return -x*sin(t) + sin(t)

for i in xrange(len(t1) - 1):
	x1[i + 1] = x1[i] + h1*f(t1[i], x1[i])
for i in xrange(len(t2) - 1):
	x2[i + 1] = x2[i] + h2*f(t2[i], x2[i])
for i in xrange(len(t3) - 1):
	x3[i + 1] = x3[i] + h3*f(t3[i], x3[i])
for i in xrange(len(t4) - 1):
	x4[i + 1] = x4[i] + h4*f(t4[i], x4[i])

from matplotlib.pyplot import plot, show, legend

plot(t1, x1, t2, x2, t3, x3, t4, x4)
legend(["1", "2", "5", "10"])
show()