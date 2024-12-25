from Oblig2_4_4c import RK4
from sys import stdout
from matplotlib.pylab import plot, xlabel, ylabel, title, legend, savefig, show, e, pi, sqrt, cos, linspace, zeros

b 		= 0.04	#kg/s
k 		= 10	#N/m
m 		= 0.1	#kg
F		= 0.1	#N
f_0	= 1 / ( 2 * pi ) * sqrt( k / m - b**2 / (2 * m**2) )
z_0		= 0.1	#m
v_0		= 0		#m/s
T		= 50	#s
N		= 1e4

f_ = linspace(f_0 / 2., f_0 * 1.5, 400)
E = zeros(len(f_))

for i in xrange(len(f_)):

	def f(v_, x_, t_):
		return - b / m * v_ - k / m * x_ + F / m * cos( 2 * pi * f_[i] * t_ )
	
	t1, z1, v1 = RK4(f, z_0, v_0, T, N)
	
	E[i] = 0.5 * m * v1[-1000: -1].max()**2

plot(f_, E)
xlabel("Frekvens")
ylabel("Energi")
title("Energi som funksjon av frekvens")
legend(["$E(t) = 0.5 m v^2$"])
savefig("Figur_%04d.png" %(7))
show()