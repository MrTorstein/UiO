#coding=utf-8
from Oblig2_i import EulerCromer

EulerCromer(dt = 0.01)
EulerCromer(dt = 0.1)

def Euler(L0 = 1, theta = 30, v0 = [0, 0], m = 0.1, k = 200., dt = 0.001):
	from numpy import linspace, zeros, sin, cos, arctan, sqrt, array, pi, linalg
	from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, savefig
	g = 9.81
	theta = theta*pi/180.
	t = 10
	n = t/dt
	r = zeros((n, 2)); r[0] = (L0*sin(theta),  - L0*cos(theta))
	v = zeros((n, 2)); v[0] = v0
	G = m*g*array([0, 1])
	for i in xrange(int(n) - 1):
		F = -G - k*(linalg.norm(r[i]) - L0)*array(r[i]/linalg.norm(r[i]))
		a = F/m
		v[i + 1] = v[i] + a*dt
		r[i + 1] = r[i] + v[i]*dt
	plot(r[:, 0], r[:, 1])
	legend(["r(t)"], loc="best")
	xlabel("x")
	ylabel("y")
	Oppgave = raw_input("Oppgave nummer: ")
	savefig("Oblig2_%s.png"%(Oppgave))
	show()

Euler(dt = 0.001)

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>Oblig2_k.py
Oppgave nummer: k1
Oppgave nummer: k2
Oppgave nummer: k3
"""