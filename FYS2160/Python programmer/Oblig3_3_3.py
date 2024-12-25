from numpy import linspace, zeros, exp, sqrt, pi, sum

T		= 300							# K
k		= 1.38 #* 10 ** (- 23)			# J / K
m		= 0.028							# kg
An		= 6.02 #* 10 ** (23)			# 1 / mol
N		= 100000						# Totalt antall steg
slutt	= 300
dv		= slutt / N

v	= linspace(0, slutt, N)				# m/s
P	= 0

for i in range( len( v ) - 1 ):
	D = ( m / ( 2 * pi * k * T * An ) ) ** ( 3 / 2 ) * 4 * pi * v[i] ** 2 * exp( - m * v[i] ** 2 / ( 2 * k * T * An ) )
	P = P + D * dv

print( "Sannsynligheten for at et nitrogen molekyl har en fart lavere enn 300 m/s er %.2f" %P )