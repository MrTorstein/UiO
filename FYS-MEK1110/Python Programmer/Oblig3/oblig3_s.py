#coding=utf-8
def txvF(m, k, v0, tslutt, u, n):
	from matplotlib.pylab import linspace, zeros_like, sign
	dt = tslutt/float(n)
	t = linspace(0, tslutt, n)
	x = zeros_like(t)
	v = zeros_like(t)
	F = zeros_like(t)
	v[0] = v0
	g = 9.81
	mys = 0.6
	myd = 0.3
	
	for i in xrange(len(t) - 1):

		if -1E-5 <= v[i] <= 1E-5:
			
			if abs(F[i]) <= mys*m*g:
				a_ = 0

			else:
				a_ = k/m*(-x[i] + u*t[i]) - myd*g*sign(v[i])

		else:
			a_ = k/m*(-x[i] + u*t[i]) - myd*g*sign(v[i])

		v[i + 1] = v[i] + a_*dt
		x[i + 1] = x[i] + v[i + 1]*dt
		F[i + 1] = k*(- x[i + 1] + u*t[i + 1])
	return t, x, v, F

if __name__ == "__main__":
	from matplotlib.pylab import plot, show, figure, legend, title, xlabel, ylabel, savefig
	figure()
	t, x, v, F = txvF(m = 0.1, k = 100., v0 = 0, tslutt = 2., u = 0.1, n = 50000)
	plot(t, x)
	legend(["x(t)"], loc = "best")
	xlabel("t")
	ylabel("x")
	title("Bevegelsen til kloss med masse 0.1kg")
	savefig("oblig3_s1.png")

	figure()
	t_, x_, v_, F_ = txvF(m = 1.0, k = 100., v0 = 0, tslutt = 2., u = 0.1, n = 50000)
	plot(t_, x_)
	legend(["x(t)"], loc = "best")
	xlabel("t")
	ylabel("x")
	title("Bevegelsen til kloss med masse 1.0kg")
	savefig("oblig3_s2.png")

	show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig3_s.py

"""