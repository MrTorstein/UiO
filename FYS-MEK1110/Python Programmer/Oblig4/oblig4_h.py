#coding=utf-8
def xt(v_0, x_0):
	U0 = 150.
	m = 23.
	x0 = 2
	alf = 39.48
	n = 1000
	t0 = 5.
	dt = t0/n
	from numpy import zeros, linspace, sign
	t = linspace(0, t0, n)
	a = zeros(n)
	v = zeros(n); v[0] = v_0
	x = zeros(n); x[0] = x_0
	for i in xrange(len(t) - 1):
		if abs(x[i]) >= x0:
			a[i] = 0
		else:
			a[i] = - U0/(m*x0)*sign(x[i]) - (alf*v[i])/m
		v[i + 1] = v[i] + a[i]*dt
		x[i + 1] = x[i] + v[i]*dt
	return x, t

if __name__ == "__main__":
	from numpy import linspace
	print "Dette programmet var ikke meningen aa kjore selv tror jeg. Det gjor hvertfall ikke noe."
	print len(linspace(8, 10, 10))
"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig4_h.py
Dette programmet var ikke meningen aa kjore selv tror jeg. Det gjor hvertfall ikke noe.

"""