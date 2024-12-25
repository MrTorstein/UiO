from numpy import linspace, zeros
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, savefig, show


def RK4(funk, tstart, tslutt, xstart, vstart, dt):
	"""
	Funksjon som tar imot initialbetingelser og en funksjon tiden skal evalueres over.
	Variabler:
		- funk		= Funksjon laget av difflikning
		- tstart	= Starttid for analysen
		- tslutt	= Sluttid for analysen
		- xstart	= Initialbetingelsen for posisjon
		- vstart	= Initialbetingelsen for hastighet
		- dt		= Endring i tid per steg.
	"""
	N	= (tslutt - tstart) / dt
	
	t	= linspace(tstart, tslutt, N)
	x	= zeros(len(t)); x[0] = xstart
	v	= zeros(len(t)); v[0] = vstart
	
	for i in xrange(len(t) - 1):
		a1	= funk(t[i], x[i], v[i])
		v1	= v[i]
		
		xh1	= x[i] + v1 * dt / 2.
		vh1 = v[i] + a1 * dt / 2.
		
		a2	= funk(t[i] + dt / 2., xh1, vh1)
		v2	= vh1
		
		xh2	= x[i] + v2 * dt / 2.
		vh2	= v[i] + a2 * dt / 2.
		
		a3	= funk(t[i] + dt / 2., xh2, vh2)
		v3	= vh2
		
		xe	= x[i] + v3 * dt
		ve	= v[i] + a3 * dt
		
		a4	= funk(t[i] + dt, xe, ve)
		v4	= ve
		
		am	= 1 / 6. * (a1 + 2 * a2 + 2 * a3 + a4)
		vm	= 1 / 6. * (v1 + 2 * v2 + 2 * v3 + v4)
		
		x[i + 1] = x[i] + vm * dt
		v[i + 1] = v[i] + am * dt
	
	return t, x, v


def diffFunk(t, x, v):
	"""
	Funksjon som gir x''(t) fra difflikninga
	m * x''(t) + k * x(t) = 0
	"""
	
	m	= 0.5	# [kg]
	k	= 1		# [N/m]
	
	a	= - k / m * x
	
	return a


if __name__ == "__main__":
	dt		= 1e-2	# [sek]
	tstart	= 0		# [sek]
	tslutt	= 20	# [sek]
	xstart	= 1		# [m]
	vstart	= 0		# [m/s]
	
	t, x, v = RK4(diffFunk, tstart, tslutt, xstart, vstart, dt)
	
	figure()
	plot(t, x)
	title("Plott av en harmonsik ocillator")
	legend(["x(t)"], loc = "upper right")
	xlabel("t [sek]")
	ylabel("x [m]")
	savefig("Prosjekt_1a.png")
	
	figure()
	plot(x, v)
	title("Plott av faserommet til funksjonen over")
	legend(["v(x)"])
	xlabel("x [m]")
	ylabel("v [m/s]")
	savefig("Prosjekt_1b.png")
	
	show()