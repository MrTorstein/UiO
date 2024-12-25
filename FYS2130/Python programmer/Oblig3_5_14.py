from numpy import linspace, sin, pi, fft, exp
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, savefig, show

Sf		= 512
frek	= 13
N		= Sf
t		= linspace(0, 1, N)
x		= 2 * sin(2 * pi * frek * t)

f	= Sf / 2 * linspace(0, 1, N / 2)
X	= fft.fft(x, N) / N

figure()
plot(t, x, "b")
xlabel("tid [s]")
ylabel("utslag [relevant enhet]")
savefig("3Figur_17.png")

figure()
plot(f, abs(X[0 : N / 2]), "r")
xlabel("Frekvens [Hz]")
ylabel("Fourierkoeffisient |X(f)|")
savefig("3Figur_18.png")

show()