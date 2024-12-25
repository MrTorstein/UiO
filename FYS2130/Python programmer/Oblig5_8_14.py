from numpy import linspace, zeros, exp, sin, cos, pi
from matplotlib.pyplot import ion, figure, draw, pause, ioff, show, axis

x0	= - 20
x1	= 20
dx	= 0.1
sig	= 2
v	= 0.5
dt	= 0.1
N	= (abs(x0) + abs(x1)) / dx + 1
x	= linspace(x0, x1, N)

def draw_wave(T):
	u = exp(- (x / (2 * sig)) ** 2)
	du_dt = (v / (2 * sig ** 2)) * x * u			# Tidsderiverte av utslaget
	u0 = u - du_dt * dt
	fac = (dt * v / dx) ** 2

	ion()
	fig = figure()
	#axis([-20, 20, -1.2, 1.2])
	ax = fig.add_subplot(111)
	line, = ax.plot(x, u)
	line2, = ax.plot(x, u)
	draw()

	for i in xrange(T):
		un = zeros(len(u))
		un[1 : - 2] = (2 * (1 - fac)) * u[1 : - 2] - u0[1 : - 2] + fac*(u[2 : - 1] + u[0 : - 3])
		un[0] = (2 * (1 - fac)) * u[0] - u0[0] + fac * u[1]
		un[- 1] = (2 * (1 - fac))*u[- 1] - u0[- 1] + fac*u[- 2]

		line.set_ydata(un)
		draw()
		pause(0.001)

		u0 = u
		u = un

	ioff()
	show()

draw_wave(1000)