from numpy import linspace, zeros, exp, sqrt, pi, sum, inf
from scipy.integrate import quad

T		= 1000							# K
k		= 1.38 #* 10 ** (- 23)			# J / K
mH		= 0.002							# kg
mHe		= 0.004							# kg
An		= 6.02 #* 10 ** (23)			# 1 / mol
start	= 11000							# m / s
slutt	= inf

def DH(v):
	return ( mH / ( 2 * pi * k * T * An ) ) ** ( 3 / 2 ) * 4 * pi * v ** 2 * exp( - mH * v ** 2 / ( 2 * k * T * An ) )

def DHe(v):
	return ( mHe / ( 2 * pi * k * T * An ) ) ** ( 3 / 2 ) * 4 * pi * v ** 2 * exp( - mHe * v ** 2 / ( 2 * k * T * An ) )

PH = quad(DH, start, slutt)
PHe = quad(DHe, start, slutt)

print( "Sannsynligheten for et hydrogen molekyl med en fart høyere enn 11 km/s, ved T = 1 kK er %.2e" %PH[0] )
print( "Sannsynligheten for et helium molekyl med en fart høyere enn 11 km/s, ved T = 1 kK er %.2e" %PHe[0] )