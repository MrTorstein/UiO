#coding=utf-8
from matplotlib.pylab import linspace, zeros, zeros_like, array, sqrt, arcsinh, cos, sin, plot, show, legend, savefig
rho = 1/2.
x = linspace(-2, 2, 1001)

s = lambda x: array((x, x**2))

XY = lambda x: array(())
def XY(x):
	a = 
	b = x**2 + rho/sqrt(1 + 4*x**2)
	return array((a, b))

def r(x):
	def m(x):
		sig = lambda x: 1/4.*(2*x*sqrt(1 + 4*x**2) + arcsinh(2*x))
		a = array((rho*(2*x*cos(sig(x)/rho) - sin(sig(x)/rho))/sqrt(1 + 4*x**2)))
		b = array((rho*(- 2*x*sin(sig(x)/rho) - cos(sig(x)/rho))/sqrt(1 + 4*x**2)))
		return array((a, b))
	return m(x) + XY(x)

plot(s(x)[0, :], s(x)[1, :])
plot(XY(x)[0, :], XY(x)[1, :])
plot(r(x)[0, :], r(x)[1, :])
legend(["s(x)", "(X(x), Y(x))", "r(x)"], loc = "best")
savefig("Oblig1_2f.png")
show()

"""
C:\Users\Torstein\Documents\UiO\Mat1110\Python programmer>Oblig1_2f.py

"""