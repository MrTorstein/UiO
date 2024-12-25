#coding=utf-8
from matplotlib.pylab import linspace, array, sqrt, pi, plot, show, scatter, polyfit, polyval
def B(x):
	mu_0 = 1.25663706*10**-6
	a = 1.75*10**-2
	t = 1*10**-2
	Js = 861818
	return mu_0/2.*Js*((x + t)/(sqrt((x + t)**2 + a**2)) - x/sqrt(x**2 + a**2))

x = linspace(0, 10**-2, 100)
y = B(x)

y2 = array([301, 57, 16.1, 6.1, 3.8, 2.6])*10**-3 - 1.1*10**-3
x2 = linspace(0, 10**-2, 6)
p = polyfit(x, y, 5)

y1 = polyval(p, x)
plot(x, y1)
scatter(x2, y2)
show()

a = 4*10**-2
t = 27**-2
x_ = 0.04
print "JS =", (2*0.8*10**-3)/(1.25663706*10**-6*(((x_ + t)/(sqrt((x_ + t)**2 + a**2)) - x_/sqrt(x_**2 + a**2))))
#1.08299
#107736.737489