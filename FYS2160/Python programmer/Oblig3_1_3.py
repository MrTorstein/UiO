from numpy import linspace, log as ln, array, zeros, arange, divide
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, axis, savefig, show

Pc	= 33.6										# [atm]
Vc	= 0.089										# [l / mol]
Tc	= 126										# [K]

T	= array([30, 77, 100, 110, 115, 120, 125])	# [K]
t	= T / Tc
v	= linspace(0.03 / Vc, 0.178 / Vc, 1E4)
p	= zeros( ( len( t ), len( v ) ) )

for i in range( len( t ) ):
	p[i]	= 8 * t[i] / ( 3 * v - 1 ) - 3 / v ** 2

leg = []

figure()
for i in range( len( t ) ):
	plot(v, p[i])
	leg.append("t = %.3f"%t[i])

title( "Plot av P - V isotermen. Dimensjonsløst uttrykk" )
legend( leg )
xlabel( "V / Vc" )
ylabel( "P / pc" )
axis( [0, 2, -10, 50] )
savefig("Figur01.png")

show()