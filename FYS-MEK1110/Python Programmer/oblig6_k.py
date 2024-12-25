#coding=utf-8
from numpy import linspace, zeros_like, array, pi, linalg, cos, sin, cross
import matplotlib.pyplot as plt

m = 0.1
k = 20.
b = 0.2
v0 = 1.
dt = 0.001

t = linspace(0, 5, 5/dt)
SM = array([zeros_like(t), zeros_like(t)]); SM[1, :] = b/2.
xA = zeros_like(SM); xA[1, 0] = b
xB = zeros_like(SM)
vA = zeros_like(SM); vA[0, 0] = v0
vB = zeros_like(SM)
aA = zeros_like(SM)
aB = zeros_like(SM)


for i in xrange(len(t) - 1):
	SM[0, i + 1] = SM[0, i] + 1./2.*v0*dt
	AB = linalg.norm(xB[:, i] - xA[:, i])
	aA[:, i] = k*(AB - b)*(xB[:, i] - xA[:, i])/(AB*m)
	aB[:, i] = -k*(AB - b)*(xB[:, i] - xA[:, i])/(AB*m)
	vA[:, i + 1] = vA[:, i] + aA[:, i]*dt
	vB[:, i + 1] = vB[:, i] + aB[:, i]*dt
	xA[:, i + 1] = xA[:, i] + vA[:, i + 1]*dt
	xB[:, i + 1] = xB[:, i] + vB[:, i + 1]*dt

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig6_k.py

"""
