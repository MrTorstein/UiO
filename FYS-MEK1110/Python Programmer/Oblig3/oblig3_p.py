#coding=utf-8
def neste(x0, v0, k, m, dt, t, u):
	a = k/m*(-x0)
	v = v0 + a*dt
	x = x0 + v*dt
	return x, v

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig3_p.py

"""