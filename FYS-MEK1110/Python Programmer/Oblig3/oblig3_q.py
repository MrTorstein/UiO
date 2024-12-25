#coding=utf-8
if __name__ == "__main__":
	from oblig3_p import neste

def vx(m, k, v0, tslutt, u, n, neste):
	from matplotlib.pylab import linspace, zeros_like
	dt = tslutt/n
	t = linspace(0, tslutt, n)
	x = zeros_like(t)
	v = zeros_like(t)
	v[0] = v0
	for i in xrange(len(t) - 1):
		x_, v_ = neste(x[i], v[i], k, m, dt, t[i], u)
		x[i + 1] = x_
		v[i + 1] = v_
	return t, x, v

if __name__ == "__main__":
	from matplotlib.pylab import plot, show, legend, xlabel, ylabel, sin, subplot, figure, savefig, sqrt
	t, x, v = vx(m = 0.1, k = 100., v0 = 0.1, tslutt = 2., n = 3000., u = 0, neste = neste)

	figure()
	subplot(211)
	omg = lambda k, m : sqrt(k/m)
	fas = 0.1/omg(100, 0.1)*sin(omg(100, 0.1)*t)
	plot(t, fas)
	xlabel("t")
	ylabel("x")
	subplot(212)
	plot(t, x)
	xlabel("t")
	ylabel("x")
	savefig("oblig3_q.png")
	show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig3_q.py

"""