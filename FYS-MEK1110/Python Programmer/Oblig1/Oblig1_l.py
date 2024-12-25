print "Lope med vind = 1m/s:"
from Oblig1_k import FC, FV		#henter FC og FV fra forrige oppgave
from Oblig1_e import Euler	#henter Euler fra forste oppgave

def F(t):
	return 400		#N

def D(t, v):
	from math import exp
	A = 0.45		#m^2
	tc = 0.67		#s
	rho = 1.293		#kg/m^3
	CD = 1.2
	w = 1			#m/s
	At = A*(1 - 0.25*exp(- (t/tc)**2))			#m^2
	return 0.5*At*rho*CD*(v - w)**2				#N

def A(t, v):
	m = 80.			#kg
	return (F(t) + FC(t) - FV(v)- D(t, v))/m	#m/s^2


t, s, v, a = Euler(A)
print "Tiden loperen bruker paa 100m er %.2fs" %t[-1]

print "\n"	#for aa faa et ekstra linjemellomrom
print "Lope mot vind = 1m/s:"

def D(t, v):
	from math import exp
	A = 0.45		#m^2
	tc = 0.67		#s
	rho = 1.293		#kg/m^3
	w = -1			#m/s
	At = A*(1 - 0.25*exp(- (t/tc)**2))			#m^2
	return 0.5*At*rho*1.2*(v - w)**2			#N


t, s, v, a = Euler(A)
print "Tiden loperen bruker paa 100m er %.2fs" %t[-1]

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig1_l.py
Lope med vind = 1m/s:
Tiden loperen bruker paa 100m er 9.21s


Lope mot vind = 1m/s:
Tiden loperen bruker paa 100m er 9.43s

"""