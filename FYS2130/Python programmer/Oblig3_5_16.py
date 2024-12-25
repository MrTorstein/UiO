from numpy import linspace, sin, pi, fft, exp
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, savefig, show

Sf		= 2**14						# Samplingsfrekvens [Hz]
frek	= 16						# Frekvens til signalet [Hz]
N		= 2**14						# Antall punkter
t		= linspace(0, 1, N)			# Tidsarray
x		= (-1)**(t // (1. / frek))	# Firkantsignalet

f	= Sf / 2 * linspace(0, 1, N / 2)
X	= fft.fft(x, N) / N

figure()
plot(t, x, "b")
xlabel("tid [s]")
ylabel("utslag [relevant enhet]")
savefig("3Figur_21.png")

figure()
plot(f, abs(X[0 : N / 2]), "r")
xlabel("Frekvens [Hz]")
ylabel("Fourierkoeffisient |X(f)|")
savefig("3Figur_22.png")

show()