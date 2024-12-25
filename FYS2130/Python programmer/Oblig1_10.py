#coding=utf-8
from matplotlib.pylab import linspace, zeros_like, plot, show, title, xlabel, ylabel, savefig

m = 1.							#kg
g = 9.81						#m/s^2
k = 1e4							#kg/s^2

t = linspace(0, 10, 1e5)		#s
v = zeros_like(t)
x = zeros_like(t); x[0] = 10	#m

dt = (t[1] - t[0])/2.
tol = 1e-5
for i in xrange(len(t) - 1):
	#print x[i]
	if x[i] <= tol:
		a = - k*x[i]/m
	else:
		a = - g
	v[i + 1] = v[i] + a*dt
	x[i + 1] = x[i] + v[i + 1]*dt

plot(x, v)
title("Faseromplott av ocillerende sprettball")
xlabel("x")
ylabel("v")
savefig("Faserom_2.png")
show()