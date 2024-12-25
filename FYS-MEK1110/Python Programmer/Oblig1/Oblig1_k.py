from numpy import exp, array
from Oblig1_i import A as Afunk	#henter A fra oppgave i, men kaller den Afunk for ikke aa forveksle med variabelen A
from Oblig1_e import Euler	#Euler fra forste oppgave
from matplotlib.pyplot import plot, legend, xlabel, show

def F(t):
	from numpy import zeros_like
	F = zeros_like(t)
	for i in xrange(len(t)):
		F[i] = 400	#N
	return F		#array med bare 400 i seg

def D(t):
	A = 0.45		#m^2
	tc = 0.67		#s
	rho = 1.293		#kg/m^3
	CD = 1.2
	At = A*(1 - 0.25*exp(- (t/tc)**2))	#m^2
	return 0.5*At*rho*CD*v**2			#N

def FC(t):
	tc = 0.67		#s
	fc = 488		#N
	return fc*exp(- (t/tc)**2)			#N

def FV(v):
	fv = 25.8		#sN/m
	return fv*v							#N

if __name__ == "__main__":		#uten betydning se tidligere oppgaver
	t, s, v, a = Euler(Afunk)
	t = array(t); v = array(v)	#gjor om lister til arrays
	
	
	plot(t, F(t), t, FC(t), t, FV(v), t, D(t))
	legend(["F(t)", "FC(t)", "FV(t)", "D(t)"], loc = "best")
	xlabel("t")
	show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>Oblig1_k.py

"""