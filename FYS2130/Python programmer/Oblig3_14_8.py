#a)
from numpy import linspace, logspace, array, meshgrid, zeros_like, sin, pi, fft, exp, cos, log, conj, absolute, floor
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, savefig, show, contourf, colorbar, subplot

Sf		= 1e4					# Samplingsfrekvens [Hz]
N		= 8192					# Antall punkter
T		= N / Sf				# Sluttid for signalet
t		= linspace(0, T, N)		# Tidsarray

f1		= 1000					# Frekvens til sinussignalet [Hz]
f2		= 1600					# Frekvens til kosinussignalet [Hz]
c1		= 1.					# Amplituden til sinussignalet [rel. enhet]
c2		= 1.7					# Amplituden til kosinussignalet [rel. enhet]

x		= c1 * sin(2 * pi * f1 * t) + c2 * sin(2 * pi * f2 * t)	# Signalet

figure()
plot(t[0 : N / 40], x[0 : N / 40], "b")	# Plotter den foorste 40ende delen av boolgen
title("To Sinus signal")
xlabel("Tid [s]")
ylabel("Utslag [relevant enhet]")
savefig("3Figur_23.png")
show()

#b)
f	= Sf * linspace(0, 1, N)
X	= fft.fft(x, N) / N

figure()
plot(f[0 : N / 2], abs(X[0 : N / 2]))
title("Sinus og kosinus etter Fouriertransformasjon")
xlabel("Frekvens [Hz]")
ylabel("Fourierkoeffisienten |X(f)|")
savefig("3Figur_24.png")
show()

#c)
def Wavelet(X, K, freks, sluttid, Fignr = 0):
	"""
	X		= Funksjon av tid som skal analyseres
	K		= Konstant som regulerer bredde paa wavelett
	freks	= Liste med to verdier, start og slutt paa analysefrekvensintervallet
	sluttid	= Tiden jeg skal plotte signalet til
	Fignr	= Nummernavn paa figuren som lagres
	"""
	M		= floor(log(freks[1]  / freks[0]) / log(1 + (1. / (8 * K)))) + 1	# Antall steg for analysefrekvensen
	omg_a	= logspace(log(freks[0]) / log(10), log(freks[1]) / log(10), M)		# Analysefrekvenser
	
	"""
	F(psi)(omg)	= 2 * (exp(- (K * (omg - omg_a)/omg_a)**2) - exp(- K**2)exp(- (K * omg / omg_a)**2))
	"""
	
	def PSI(omg, omg_a):
		"""
		omg		= Array av frekvenser brukt i Fouriertransformasjonen av signalet
		omg_a	= Array med alle frekvensene jeg vil analysere over (y-aksen i plottet)
		"""
		return conj(2 * (exp(- (K * (omg - omg_a)/omg_a)**2) - exp(- K**2) * exp(- (K * omg / omg_a)**2)))
	
	M = []
	for i in xrange(len(omg_a)):
		M.append(fft.ifft(X * PSI(f, omg_a[i])))
	
	M = array(M[0 : sluttid][0 : sluttid])
	T_, O_ = meshgrid(t[0 : sluttid], omg_a[0 : sluttid])
	
	
	figure()
	ax = subplot()
	ax.set_yscale("log")
	A = contourf(T_, O_, absolute(M[:, 0 : sluttid]))
	colorbar(A)
	title("Konturplott av Wavelet analysen med K = %d"%K)
	xlabel("Tid [s]")
	ylabel("Frekvens [Hz]")
	savefig("3Figur_%02d.png"%Fignr)

Wavelet(X, 24, [800, 2000], N, 25)
Wavelet(X, 200, [800, 2000], N, 26)
show()

#d)
t1		= 0.15					# Konstant [s]
t2		= 0.5					# Konstant [s]
siggy1	= 0.01					# Konstant [s]
siggy2	= 0.1					# Konstant [s]

f_t = c1 * sin(2 * pi * f1 * t) * exp(- ((t - t1) / siggy1)**2) + c2 * sin(2 * pi * f2 * t) * exp(- ((t - t2) / siggy2)**2)

figure()
plot(t, f_t, "b")	# Plotter boolgen
title("Sinus og eksponensial signal")
xlabel("Tid [s]")
ylabel("Utslag [relevant enhet]")
savefig("3Figur_27.png")


F_T = fft.fft(f_t, N) / N

figure()
plot(f[0 : N / 2], abs(F_T[0 : N / 2]))
title("Sinus og eksponensial etter Fouriertransformasjon")
xlabel("Frekvens [Hz]")
ylabel("Fourierkoeffisienten |X(f)|")
savefig("3Figur_28.png")
show()


Wavelet(F_T, 24, [800, 2000], N, 29)
Wavelet(F_T, 100, [800, 2000], N, 30)
Wavelet(F_T, 8, [800, 2000], N, 31)
Wavelet(F_T, 200, [800, 2000], N, 32)
show()