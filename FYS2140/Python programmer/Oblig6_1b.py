from numpy import linspace, exp, sqrt, pi
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig

m	= 1
omg	= 1
h_	= 1.0545718 * 10 ** ( - 3 )

fratil = 0.2
x		= linspace( - fratil, fratil, 1e5 )
psi_0	= lambda x: ( ( m * omg ) / ( pi * h_ ) ) ** ( 1 / 4. ) * exp( - ( m * omg ) / ( 2 * h_ ) * x ** 2 )
psi_1	= lambda x: psi_0( x ) * sqrt( ( 2 * m * omg ) / ( h_ ) ) * x
psi_2	= lambda x: sqrt( 2 ) * psi_0( x ) * ( ( 2 * m * omg * x ** 2 ) / ( h_ ) - 1 )

plot(x, psi_0( x ), x, psi_1( x ), x, psi_2( x ))
legend(["$ \psi_0(x) $", "$ \psi_1(x) $", "$ \psi_2(x) $"])
title("Plott av de tre foorste HO-ene")
xlabel("Posisjon")
ylabel("Sannsynlighet")
savefig("Figur_04.png")
show()
