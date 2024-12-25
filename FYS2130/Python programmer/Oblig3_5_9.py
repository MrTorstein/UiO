# Enkelt eksempelprogram for aa vise hvordan fouriertransformasjon
# kan gjennomfoores i praksis i python. Eksempelet er en modifikasjon
# av et eksempelprogram paa hjemmesidene til Matlab og endret av meg
# til pythonkode. Kan ikke skrive de tre Norske bokstavene, derfor har
# jeg endret dem til ae, oo og aa.

from numpy import linspace, zeros, pi, cos, fft
from matplotlib.pyplot import plot, title, figure, xlabel, ylabel, savefig, show

Fs	= 1000									# Samplingsfrekvensen
dt	= 1. / Fs								# Tid mellom hver sampling
N	= 1024									# Antall samplinger
t	= linspace(0, N - 1, N) * dt			# Tidsvektor

# Lager her et kunstig signal som er et enkelt kosinus signal.
# Legger inn frekvensen som en variabel slik at jeg enkelt kan variere denne.
frek	= 700.								# Frekvens i Hz
x		= 0.8 * cos(2 * pi * frek * t)		# Signal, en enkel kosinus

figure()									# Lager figur
plot(Fs * t, x)								# Plott av signalet i tidsbildet
title("Opprinnelig signal (tidsbilde)")
xlabel("tid (millisekunder)")
savefig("3Figur_00.png")					# Lagrer plottet

X = fft.fft(x, N) / N						# Fouriertransformasjon

frekv = (Fs / 2.) * linspace(0, 1, N / 2)	# Frekvensvektor for plott

# Plotter bare lengden paa frekvenskomponentene i frekvensspekteret.
# Velger aa bare ta med frekvenser opp til halve samplingsfrekvensen.
figure()									# Lager ny figure
plot(frekv, 2 * abs(X[0 : N / 2]))			# Plotter halvparten av Fourierspekteret
title("Absolutt-verdier av frekvensspekteret")
xlabel("Frekvens (Hz)")
ylabel("|X(frekv)|")
savefig("3Figur_01.png")					# Lagrer plottet
show()