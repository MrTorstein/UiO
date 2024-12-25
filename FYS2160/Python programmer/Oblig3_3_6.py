from numpy import linspace, zeros, exp, sqrt, pi, sum, inf
from scipy.integrate import quad

T		= 1000							# K
k		= 1.38 #* 10 ** (- 23)			# J / K
m		= 0.028							# kg
An		= 6.02 #* 10 ** (23)			# 1 / mol
start	= 2400							# m / s
slutt	= inf

def D(v):
	return ( m / ( 2 * pi * k * T * An ) ) ** ( 3 / 2 ) * 4 * pi * v ** 2 * exp( - m * v ** 2 / ( 2 * k * T * An ) )

P = quad(D, start, slutt)

print( "Sannsynligheten for et nitrogen molekyl med en fart høyere enn 2.4 km/s, ved T = 1 kK er %.2e" %P[0] )