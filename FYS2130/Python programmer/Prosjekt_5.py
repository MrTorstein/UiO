from Prosjekt_1 import RK4
from numpy import cos, sqrt, sin
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, savefig, show

def diffFunk(t, x, v):
	"""
	Funksjon som gir x''(t) fra difflikninga
	m * x''(t) + b * x'(t) + k * x(t) = F_D * cos(omg_D * t)
	"""
		
	m		= 0.5				# [kg]
	b		= 0.1				# [kg/s]
	k		= 1					# [N/m]
	F		= 0.7				# [N]
	omg0	= sqrt(k / m)		# Svingefrekvensen til fri harmonsik ocillator
	omgD	= 13 / (8) * omg0	# Drivefrekvensen
	
	a	= - k / m * x - b / m * x + F * cos(omgD * t)
	
	return a


if __name__ == "__main__":
	dt		= 1e-2	# [sek]
	tstart	= 0		# [sek]
	tslutt	= 100	# [sek]
	xstart	= 2		# [m]
	vstart	= 0		# [m/s]
	
	
	t, x, v = RK4(diffFunk, tstart, tslutt, xstart, vstart, dt)
	
	figure()
	plot(t, x)
	title("Plott av en harmonsik ocillator med paatrykt kraft og dempning")
	legend(["x(t)"])
	xlabel("t [sek]")
	ylabel("x [m]")
	savefig("Prosjekt_5a.png")
	
	figure()
	plot(x, v)
	title("Plott av faserommet til funksjonen over")
	legend(["v(x)"])
	xlabel("x [m]")
	ylabel("v [m/s]")
	savefig("Prosjekt_5b.png")
	
	
	show()