from numpy import linspace, zeros, pi
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, savefig, show


def RK4(funk, tstart, tslutt, xstart, vstart, dt, m):
	"""
	Funksjon som tar imot initialbetingelser og modelerer en draape over tid.
	Variabler:
		- funk		= Funksjon laget av difflikning
		- tstart	= Starttid for analysen
		- tslutt	= Sluttid for analysen
		- xstart	= Initialbetingelsen for posisjon
		- vstart	= Initialbetingelsen for hastighet
		- dt		= Endring i tid per steg.
		- m			= Utgangsmassen til draapen
	"""
	
	x_c	= 0.0025	# [m]
	B	= 50		# [s/m]
	psi	= 0.00055	# [kg/s]
	rho	= 1000		# [kg/m^3]

	N	= (tslutt - tstart) / dt
	
	t	= linspace(tstart, tslutt, N)
	x	= zeros(len(t)); x[0] = xstart
	v	= zeros(len(t)); v[0] = vstart
	
	for i in xrange(len(t) - 1):
		a1	= funk(t[i], x[i], v[i], m)
		v1	= v[i]
		
		xh1	= x[i] + v1 * dt / 2.
		vh1	= v[i] + a1 * dt / 2.
		
		a2	= funk(t[i] + dt / 2., xh1, vh1, m)
		v2	= vh1
		
		xh2	= x[i] + v2 * dt / 2.
		vh2	= v[i] + a2 * dt / 2.
		
		a3	= funk(t[i] + dt / 2., xh2, vh2, m)
		v3	= vh2
		
		xe	= x[i] + v3 * dt
		ve	= v[i] + a3 * dt
		
		a4	= funk(t[i] + dt, xe, ve, m)
		v4	= ve
		
		am	= 1 / 6. * (a1 + 2 * a2 + 2 * a3 + a4)
		vm	= 1 / 6. * (v1 + 2 * v2 + 2 * v3 + v4)
		
		x[i + 1] = x[i] + vm * dt
		v[i + 1] = v[i] + am * dt
		
		m	= m + psi * dt
	
		if x[i + 1] > x_c:
			dm	= B * m * v[i]
			dx	= ((3 * dm ** 4) / (4 * pi * rho * m ** 3)) ** (1 / 3.)
			m = m - dm
			
			x[i + 1] = x[i + 1] - dx
	
	return t, x, v


def diffFunk(t, x, v, m):
	"""
	Funksjon som gir x''(t) fra difflikninga
	m(t) * x''(t) + (b + psi) * x'(t) + k * x(t) = m(t) * g
	"""
	
	psi	= 0.00055	# [kg/s]
	b	= 0.001		# [kg/s]
	k	= 0.475		# [N/m]
	g	= 9.81		# [m/s^2]
	
	
	a	= m * g - k / m * x - (b + psi) / m * v
	
	return a


if __name__ == "__main__":
	dt		= 1e-4		# [sek]
	tstart	= 0			# [sek]
	tslutt	= 21		# [sek]
	xstart	= 0.001		# [m]
	vstart	= 0.001		# [m/s]
	mstart	= 0.00001	# [kg]
	
	t, x, v = RK4(diffFunk, tstart, tslutt, xstart, vstart, dt, mstart)
	
	figure()
	plot(t, x * 1e3)
	title("Plott av en model av dryppende kran")
	legend(["x(t)"])
	xlabel("t [sek]")
	ylabel("x [m] * 10^-3")
	savefig("Prosjekt_7a.png")
	
	figure()
	plot(x, v)
	title("Plott av faserommet til modellen av dryppende kran")
	legend(["v(x)"])
	xlabel("x [m]")
	ylabel("v [m/s]")
	savefig("Prosjekt_7b.png")
	
	show()