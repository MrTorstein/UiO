""" Program for running the Energy class and generating plots """
from Prosjekt1 import Energy
from numpy import logspace, zeros, sum as Sum
from matplotlib.pyplot import figure, plot, show, legend, xlabel, ylabel, savefig, xscale, rcParams

# Øker skriftstørrelsen på plottet #
rcParams.update({"font.size": 22})

# Lager temperatur liste mellom 10^4 og 10^9 #
Temp = logspace(4, 9, 10000)
Dens = 1.62 * 10 ** 5

Res = zeros([len(Temp), 4])

# Kjører utregning og lagrer resultater #
for i in range(len(Temp)):
    A = Energy(Temp[i], Dens, "N")
    Res_tot = A.Kalk()
    Res[i]  = A.eps[0:4]

# Plotter resultatene #
figure(figsize = (8, 7))

for i in range(4):
    plot(Temp, Res[:, i] / Sum(Res, 1))

legend([r"PPI", r"PPII", r"PPIII", r"CNO"])
xscale("log")
xlabel(r"Temp [K]")
ylabel(r"$\epsilon / \epsilon_{Tot}$")
savefig("Figur_00.png")
show()
