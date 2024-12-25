from Oblig6_15_15 import HvitStoyGauss
from numpy import linspace, logspace, array, meshgrid, zeros_like, sin, pi, fft, exp, cos, log, conj, absolute, floor
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, savefig, show, contourf, colorbar, subplot

def Fourier(x, Sf, Fignr):
	"""
	x		= Signalet
	Sf		= Samplingsfrekvensen
	"""
	N	= len(x)
	f	= Sf * linspace(0, 1, N)
	X	= fft.fft(x, N) / N

	figure()
	plot(f[0 : N / 2], abs(X[0 : N / 2]))
	title("Signalet etter Fouriertransformasjon")
	xlabel("Frekvens [Hz]")
	ylabel("Fourierkoeffisienten |X(f)|")
	savefig("Figur_%04d.png"%Fignr)
	show()
	
	return f, X

def Wavelet(f, X, K, freks, sluttid, Fignr = 0):
	"""
	f		= Frekvensliste fra 0 til samplingsfrekvensen
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
	savefig("Figur_%04d.png"%Fignr)
	show()

if __name__ == "__main__":
	Sf		= 2
	senfre	= 5000
	FBf		= 3000
	N		= 4000
	t, y, q = HvitStoyGauss(Sf, N, senfre, FBf)
	
	f, Y	= Fourier(y, Sf, 9)
	Wavelet(f, Y, 24, [0.1, 3], int(t[-1]) + 1, 10)
	
	y_		= y ** 2
	f, Y	= Fourier(y, Sf, 11)
	Wavelet(f, Y, 24, [0.1, 3], int(t[-1]) + 1, 12)