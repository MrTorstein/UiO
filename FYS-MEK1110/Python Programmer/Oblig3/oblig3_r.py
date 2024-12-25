#coding=utf-8
def neste(x0, v0, k, m, dt, t, u):
	a = k/m*(-x0 + u*t)
	v = v0 + a*dt
	x = x0 + v*dt
	return x, v

if __name__ == "__main__":
	from oblig3_q import vx
	from matplotlib.pylab import plot, show, legend, xlabel, ylabel, savefig, sqrt, sin, cos
	t, x, v = vx(m = 0.1, k = 100., v0 = 0, tslutt = 2., n = 3000., u = 0.1, neste = neste)
	fas = 0.1*t - 0.1/sqrt(100./0.1)*sin(sqrt(100/0.1)*t)
	plot(t, x, t, fas)
	legend(["numerisk", "analytisk"], loc = "best")
	xlabel("t")
	ylabel("x")
	savefig("oblig3_r.png")
	show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig3_r.py

"""