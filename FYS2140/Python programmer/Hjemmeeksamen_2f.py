from numpy import *
from numpy.ma import *
from matplotlib.pyplot import *
E0p		= 3727			# Hvileenergi for proton [MeV]
hbarc	= 0.1973		# [MeV pm]
c		= 3.0e2			# Lyshastighet [pm / as]
x_1		= 7.3e-3		# [pm]
x_2		= x_1 + 10e-3	# [pm]
def Psi0( x ) :
	"""
	Initialtilstand for gaussisk boelgepakke. Deler formelen
	opp i konstant og uttrykk for aa faa det mer oversiktlig.
	"""
	x0	= 0.005	# [pm]
	a	= 0.001 # [pm]
	k	= 1380	# [ 1 / pm]
	A	= ( 1. / ( 2 * pi * a ** 2 ) ) ** 0.25			# [1 / sqrt( pm )]
	K1	= exp( - ( x - x0 ) ** 2 / ( 4. * a ** 2 ) )	# [ Enhetsloos ]
	K2	= exp( 1j * k * x )								# [ Enhetsloos ]
	return A * K1 * K2									# [1 / sqrt( pm )]

def V( x ) :
	""" En potensial barriere """
	# Bruker masking for aa sette alle verdier inni barrieren lik 1
	return 34 * ( masked_inside(x, x_1, x_2).mask )

N_x	= 800	# Antall punkter i x-retning
N_p	= 1e2	# Plotte kun hver hver N_x'te utregning ( bedre "framerate" )
dx	= 0.001	# Avstand mellom x-punkter

# Lar halvparten av x-verdiene variere paa hver sin side av null
a = -0.5 * N_x * dx
b = 0.5 * N_x * dx
x = linspace( a, b, N_x )

# Erklaerer tomme arrayer for den deriverte og Psi .
# Spesifiserer at den deriverte er kompleks
dPsidx2 = zeros( N_x ).astype( "complex64" )
Psi = Psi0( x )

# Tidsparametre
T = 0.002 # Hvor lenge simuleringen skal kjoere [ as ]
dt = 1e-7 # Avstand mellom tidssteg [ as ]

# Henter ut potensialet
V_x = V( x )

# Aktiverer interaktivmodus
ion()

# Tegner initialtilstand
figure()
title("Plot av Boolgefunksjon, avhengig av tid")
xlabel("Posisjon [pm]")
ylabel("Utslag [1 / pm]")
plot([x_1 - 1e-10, x_1, x_2, x_2 + 1e-10], [0, 200, 200, 0], "--")
axis((-0.03, 0.03, 0, 600))
line , = plot( x, abs( Psi ) ** 2 )
draw()

#Konstanter
c1 = ( 1j * hbarc * c ) / ( 2. * E0p )	# [pm^2 / as]
c2 = - ( 1j * c ) / hbarc				# [1 / as MeV]
t = 0 # Teller tid
c = 1 # Teller antall utregninger

while t < T:
	# Regner ut den deriverte
	dPsidx2 [ 1 : N_x-1] = ( Psi [ 2 : N_x ] - 2 * Psi [ 1 : N_x-1 ] + Psi [ 0 : N_x-2 ] ) / dx ** 2
	
	# Regne ut ny Psi med potensialet
	Psi = Psi + dt * ( c1 * dPsidx2 + c2 * V_x * Psi )
	
	# Sjekker om den skal plotte denne utregningen
	if c == N_p :
		line.set_ydata( abs( Psi ) ** 2 )
		draw()
		pause(1e-40)
		c = 0 # Nullstillerteller
		legend(["Potensialbarrieren", "$|\Psi(x)|^2$ ved t = %02f [as]"%(t)], loc=1)
	
	# Lager bilder foor, under og etter interagering med potensial
	if 0.00190 < t < 0.0019001:
		savefig("Figur_HE_3.png")
	elif 0.00020 < t < 0.000201:
		savefig("Figur_HE_2.png")
	elif 0.00001 < t < 0.000011:
		savefig("Figur_HE_1.png")
	
	t += dt # Legger til et tidssteg
	c += 1 # Legger til en utregning

# Av med interaktiv modus
ioff()

# Beholder vinduene
show()