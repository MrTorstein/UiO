def A(t, v):
	from math import exp
	F = 400.		#N
	m = 80.			#kg
	fv = 25.8		#sN/m
	fc = 488		#N
	A = 0.45		#m^2
	tc = 0.67		#s
	rho = 1.293		#kg/m^3
	CD = 1.2
	At = A*(1 - 0.25*exp(- (t/tc)**2))	#m^2
	D = 0.5*At*rho*CD*v**2				#N
	FC = fc*exp(- (t/tc)**2)			#N
	FV = fv*v							#N
	return (F + FC - FV - D)/m			#m/s^2

if __name__ == "__main__":
	from Oblig1_e import Euler	#bruker Euler fra forste oppgave
	from matplotlib.pyplot import plot, legend, xlabel, show, savefig
	
	t, s, v, a = Euler(A)
	plot(t, s, t, v, t, a)
	legend(["s(t)", "v(t)", "a(t)"], loc = "best")
	xlabel("t")
	savefig("Oblig1_i.png")
	show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>Oblig1_i.py

"""