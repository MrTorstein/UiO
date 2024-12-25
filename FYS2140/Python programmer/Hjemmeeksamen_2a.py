from numpy import linspace, sqrt, pi, e, imag
from matplotlib.pyplot import plot, legend, title, xlabel, ylabel, savefig, show

x_0 = 5		# [fm]
a	= 1		# [fm]
k	= 1.38	# [fm^(-1)]

x	= linspace(0, 10, 1e4)
psi	= 1 / sqrt( 2 * pi * a ** 2 ) * ( e ** ( - ( ( x - x_0 ) ** 2 ) / ( 4 * a ** 2 ) ) ) ** 2

plot( x, psi )
legend(["$|\Psi(x, 0)|^2$"])
title("Plott av absoluttverdien til boolgefunksjonen, i annen")
xlabel("Posisjon [fm]")
ylabel("Utslag [1 / fm]")
savefig("Figur_06.png")
show()