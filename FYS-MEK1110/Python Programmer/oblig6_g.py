#coding=utf-8
from numpy import linspace, zeros_like, array
import matplotlib.pyplot as plt

m = 0.1
k = 20.
b = 0.2
v0 = 1.
dt = 0.001

t = linspace(0, 5, 5/dt)
x = array([zeros_like(t), zeros_like(t)]); x[0, 0] = - b
v = zeros_like(x); v[0, 0] = v0
a = zeros_like(x)
SM = zeros_like(t); SM[0] = - b/2.

for i in xrange(len(t) - 1):
	a[:, i] = [k*(x[1, i] - x[0, i] - b)/m, -k*(x[1, i] - x[0, i] - b)/m]
	v[:, i + 1] = v[:, i] + a[:, i]*dt
	x[:, i + 1] = x[:, i] + v[:, i + 1]*dt
	SM[i + 1] = SM[i] + 1./2.*v0*dt

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig6_g.py

"""