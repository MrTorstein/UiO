#coding=utf-8
from oblig4_h import xt
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, title
from sys import exit
from numpy import linspace
s = 0
tol = 1E-10
liste = linspace(8, 10, 1000)
for i in liste:
	x, t = xt(i, -5)
	if s == 1:
		exit()

	elif (-2 + tol) < x[-1] < (2 + tol):
		pass

	elif (-2 - tol) < x[-1] < (2 - tol):
		pass

	else:
		print i - (10 - 8)/1000., "slapp unna =("
		s = 1

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig4_k.py
8.66066266266 slapp unna =(

"""