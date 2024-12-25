from numpy import array, linspace, zeros, exp, pi, sum as Sum
from cmath import sqrt as csqrt

Signal = array([1, 3, 1, 3, 1, 3, 1, 3, 1, 3])

def FT(x):
	x = array(x, dtype = "complex")
	N = len(x)
	n = linspace(0, N - 1, N, dtype = "complex")
	liste = zeros(N, dtype = "complex")
	
	for i in xrange(len(n)):
		liste[i] = 1. / N * Sum(x * exp(- csqrt(- 1) * 2 * pi * n[i] * n / N))
	return liste

X = FT(Signal)

print X[0], Sum(Signal) / len(Signal)

"""
C:\Users\Torstein\Documents\UiO\Fys2130\Python programmer>python Oblig3_5_8.py
(2+0j) 2

"""