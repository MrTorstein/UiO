from numpy import linspace, log as ln, array, zeros, arange, divide
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, axis, savefig, show

Pc	= 33.6				    			# [atm]
Vc	= 0.089	    						# [l / mol]
Tc	= 126			    				# [K]

T	= 30	    						# [K]
t	= T / Tc
v	= linspace(0.03 / Vc, 100 / Vc, 1E5)

p	= 8 * t / ( 3 * v - 1 ) - 3 / v ** 2

figure()
plot( p, v, [0, 1.E-10], [0, 1000] )
title( "Plot av P - V isotermen til Maxwell equal area" )
legend( ["t = %.3f"%t] )
ylabel( "V / Vc" )
xlabel( "P / pc" )
axis( [-0.1, 0.1, 0, 1200] )
savefig("Figur02.png")

show()

#Dette er en test