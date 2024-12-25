from numpy import random, sum as Sum
from matplotlib.pyplot import figure, plot, hist, title, xlabel, ylabel, savefig, show

M = 10**4						# Antall mikrotilstander
N = 60							# Antall forskjellige spinn
d = 10**2						# Mellomrom mellom tilstander som blir plottet

liste = random.randint(-1, 1, (M, N)) * 2 + 1	# Lager en M X N liste av 1 eller -1
summer = Sum(liste, 1)							# Summerer spinnet til mikrotilstandene

#Plotter og lager histogram
figure()
title("Plott av hver hundrede mikrotilstand sitt spinn")
xlabel("Mikrotilstand / 100")
ylabel("Spinn")
plot(range(0, d), summer[0: -1: M / d])
savefig("Figur102.png")

figure()
title("Histogram av antall mikrotilstander med ett spinn")
xlabel("Spinn")
ylabel("Antall mikrotilstander")
hist(summer)
savefig("Figur103.png")

show()