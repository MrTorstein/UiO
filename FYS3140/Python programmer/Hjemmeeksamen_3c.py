from matplotlib.pylab import linspace, zeros, sin, figure, plot, legend, title, savefig, show

N = int(1e5)							# Antall punkter
n = 20									# presisjonen til delta-sekvensen

def phi(x, n):							# Definerer phi_n
	if abs(x) >= 1. / n:
		return 0
	else:
		return n / 2.

t = linspace(-2, 2, N)					# Punktene jeg evaluerer phi_n over
Phi	= zeros(N)

a = 1									# Velger en konstant
f_i = t**2 - a**2						# Evaluerer funksjon nummer 1 over t

for i in xrange(len(t)):				# Putter dette inn i phi_n
	Phi[i] = phi(f_i[i], n)

figure()								# Plotter
plot(t, Phi)
legend(["Phi_%d(t^2 - a^2)"% n])
title("Deltasekvens av funksjonen i 3b.i)")
savefig("Figur_HE_1.png")

t = linspace(-4, 4, N)					# Velger et nytt intervall
f_ii = sin(t)							# Evaluerer funksjon 2 over ny t

for i in xrange(len(t)):				# Putter ny funksjon inn i phi_n
	Phi[i] = phi(f_ii[i], n)

figure()								# Plotter 
plot(t, Phi)
legend(["Phi_%d(sin(t))"% n])
title("Deltasekvens av funksjonen i 3b.ii)")
savefig("Figur_HE_2.png")

show()									# Viser figurer