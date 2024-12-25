#coding=utf-8
from matplotlib.pylab import linspace, zeros_like, plot, show, title, xlabel, ylabel, savefig

m = 1.							#kg
g = 9.81						#m/s^2
k = 0.5							#kg/s^2

t = linspace(0, 20, 100000)		#s
v = zeros_like(t)
x = zeros_like(t)

dt = (t[1] - t[0])/2.
for i in xrange(len(t) - 1):
	a = (- k*x[i] - m*g)/m
	v[i + 1] = v[i] + a*dt
	x[i + 1] = x[i] + v[i + 1]*dt

plot(x, v)
title("Faseromplott av ocillerende fjaer og lodd")
xlabel("x")
ylabel("v")
savefig("Faserom_1.png")
show()