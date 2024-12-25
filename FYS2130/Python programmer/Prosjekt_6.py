from Prosjekt_1 import RK4
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, savefig, show

def diffFunk(t, x, v):
	"""
	Funksjon som gir x''(t) fra difflikninga
	m(t) * x''(t) + (b + psi) * x'(t) + k * x(t) = m(t) * g
	"""
	
	psi	= 0.00055	# [kg/s]
	b	= 0.001		# [kg/s]
	k	= 0.475		# [N/m]
	g	= 9.81		# [m/s^2]
	m	= 0.00001	# [kg]
	
	m_	= m + psi * t
	
	a	= m_ * g - k / m_ * x - (b + psi) / m_ * v
	
	return a


if __name__ == "__main__":
	dt		= 1e-4	# [sek]
	tstart	= 0		# [sek]
	tslutt	= 3		# [sek]
	xstart	= 0.001	# [m]
	vstart	= 0.001	# [m/s]
	
	t, x, v = RK4(diffFunk, tstart, tslutt, xstart, vstart, dt)
	
	figure()
	plot(t, x * 1e3)
	title("Plott av en harmonsik ocillator med dempning")
	legend(["x(t)"], loc = "upper right")
	xlabel("t [sek]")
	ylabel("x [m] * 10^-3")
	savefig("Prosjekt_6a.png")
	
	figure()
	plot(x, v)
	title("Plott av faserommet til funksjonen over")
	legend(["v(x)"])
	xlabel("x [m]")
	ylabel("v [m/s]")
	savefig("Prosjekt_6b.png")
	
	show()