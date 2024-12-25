from numpy import linspace, zeros, exp, transpose, cos, sin, conj, floor, sum as Sum, pi, round
from numpy.random import rand
from numpy.fft import ifft
from cmath import sqrt as csqrt
from matplotlib.pyplot import plot, show, savefig, xlabel, ylabel

def HvitStoyGauss(Fs, N, fsenter, fullFbredde):
	"""
	Parametre: 
	Fs			: Samplingsfrekvens
	N			: Antall datapunkt
	fsenter		: senterfrekvens
	fullFbredde	: Frekvensspekteret har en gaussisk fordeling
	med senterfrekvens fsenter og full bredde (1/e) i
	frekvensspekteret lik fullFbredde..
	"""
	fsigma = fullFbredde / 2.
	y = zeros(N)
	T = N / Fs
	t = linspace(0, T * (N - 1) / N, N)
	f = linspace(0, Fs * (N - 1) / N, N)
	nsenter = floor(N * fsenter / (Fs * (N - 1) / N))
	nsigma = floor(N * fsigma / (Fs * (N - 1) / N))
	gauss = exp(- (f - fsenter) * (f - fsenter) / (fsigma * fsigma))
	ampl = rand(N)
	ampl = ampl * transpose(gauss)
	faser = rand(N)
	faser = faser * 2 * pi
	y = ampl * (cos(faser) + csqrt(- 1) * sin(faser))
	
	#Speiler nedre del rundt (Nhalv + 1) for aa faa ovre del korrekt
	Nhalv = round(N / 2)
	
	for k in xrange(1, Nhalv - 1):
		y[N - k] = conj(y[k + 1])
		
	y[Nhalv + 1] = y[Nhalv + 1].real
	y[0] = 0
	q = (ifft(y) * 200).real
	return t, y, q


def autokorr(g, M):
	N = len(g)
	C = zeros(N)
	for j in xrange(N - M):
		teller = 0
		nevner = 0
		for i in xrange(M):
			teller += g[i] * g[i + j]
			nevner += g[i] * g[i]
		
		C[j + 1] = teller / nevner
	return C

if __name__ == "__main__":
	t, y, q = HvitStoyGauss(2, 4000, 5000, 3000)
	f = autokorr(y, 2000)
	
	plot(t, f)
	xlabel("Tid [sek]")
	ylabel("Autokorrelasjon [relevant enhet]")
	savefig("Figur_0008.png")
	show()