from Prosjekt_1 import RK4
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, savefig, show

def diffFunk(t, x, v):
	"""
	Funksjon som gir x''(t) fra difflikninga
	m * x''(t) + b * x'(t) + k * x(t) = 0
	"""
		
	m	= 0.5	# [kg]
	b	= 0.1	# [kg/s]
	k	= 1		# [N/m]
	
	a	= - k / m * x - b / m * v
	
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
	title("Plott av en harmonsik ocillator med dempning")
	legend(["x(t)"], loc = "upper right")
	xlabel("t [sek]")
	ylabel("x [m]")
	savefig("Prosjekt_2a.png")
	
	figure()
	plot(x, v)
	title("Plott av faserommet til funksjonen over")
	legend(["v(x)"])
	xlabel("x [m]")
	ylabel("v [m/s]")
	savefig("Prosjekt_2b.png")
	
	show()