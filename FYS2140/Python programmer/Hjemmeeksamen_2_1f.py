from numpy import pi, sqrt, exp, linspace, zeros
from numpy.ma import masked_inside
from matplotlib.pyplot import figure, ion, ioff, draw, plot, show, legend, axis, xlabel, ylabel, title, pause

""" Variabler """
h_c		= 1.973e2	# [MeV fm]
E0p		= 3727		# [MeV]
c		= 3e5		# [fm / as]
x_0		= 5			# [fm]
a		= 1			# [fm]
k		= 1.38		# [fm^-1]
V_0		= 34		# [MeV]
A		= ( 1. / ( 2 * pi * a ** 2 ) ) ** ( 1. / 4 )
x_1		= 7.3		# [fm]
x_2		= x_1 + 10	# [fm]
T		= 1e-3		# [as]
N_x		= int(1e3)	# Antall x verdier
N_t		= int(1e4)	# Antall t verdier
N_c		= int(1e2)	# Avstand mellom utregninger som skal plottes
dx		= ( 2 * ( x_2 + 1 ) ) / N_x
dt		= T / N_t

x	= linspace(- 1 - x_2, 1 + x_2, N_x)
Psi	= zeros(N_x).astype("complex64")
a_x	= zeros(N_x).astype("complex64")

Psi0	= lambda x: A * exp( - ( x - x_0 ) ** 2 / ( 4 * a ** 2 ) ) * exp( 1j * k * x )
V		= lambda x: V_0 * masked_inside(x, x_1, x_2).mask

Psi = Psi0(x)
V_x = V(x)

"""
Schroodinger likningen:
d Psi / dt = ( h_ i / 2 m ) ( d^2 Psi / dx^2 ) - ( V i / h_ ) Psi
a(1) = v(1.5) - v(0.5) / dx = ( psi(2) - 2 * psi(1) + psi(0) ) / ( dx ** 2 )
v(1.5) = psi(2) - psi(1) / dx
"""

ion()
figure()
axis((x[0], x[-1] + 1, 0, 5))
plot([x_1 - 1e-10, x_1, x_2, x_2 + 1e-10], [0, 1, 1, 0])
linje = plot(x, abs(Psi) ** 2)
draw()

t = 0
while t < T:
	a_x[1 : N_x - 1] = ( Psi[2 : N_x] - 2 * Psi[1 : N_x - 1] + Psi[0 : N_x - 2] ) / dx ** 2
	v = ( h_c * c * 1j / ( 2. * E0p ) ) * a_x - ( V_x * 1j * c / h_c ) * Psi
	Psi = Psi + v * dt
	if c == N_c:
		linje[0].set_ydata(abs(Psi) ** 2)
		draw()
		pause(1e-10)
	
	t += dt
	
ioff()

show()