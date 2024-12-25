from numpy import array, linspace, zeros, pi, exp, sin, cos, sum as Sum, fft
from matplotlib.pyplot import figure, plot, title, legend, xlabel, ylabel, savefig, show

with open("SN_y_tot_V2.0.txt", "r") as innfil:
	tider = []
	solflekker = []
	for linje in innfil:
		koll = linje.split()
		tider.append(float(koll[0]))
		solflekker.append(float(koll[1]))
	tider = array(tider)
	solflekker = array(solflekker)

# Plotter tidsbildet
figure()
plot(tider, solflekker, "b")
title("Tidsbilde")
xlabel("Tid [Aar]")
ylabel("Antall solflekker")
savefig("3Figur_15.png")

N		= len(tider)							# Antall samplinger
Sf		= N / (tider[-1] - tider[0])			# Samplingsfrekvensen
frek	= (Sf / 2) * linspace(0, 1, N / 2)		# Frekvensarray

X		= fft.fft(solflekker, N) / N			# Fouriertransformasjonen

figure()										# Plotter frekvensbildet
plot(frek, abs(X[0 : N / 2]))
title("Frekvensbilde")
xlabel("Frekvens [$aar^{-1}$]")
ylabel("Fourierkoeffisienten (|X(f)|)")
savefig("3Figur_16.png")

show()