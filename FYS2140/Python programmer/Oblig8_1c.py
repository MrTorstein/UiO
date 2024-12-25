from numpy import linspace, sqrt
from matplotlib.pyplot import plot, legend, title, xlabel, ylabel, savefig, show

E = linspace(1, 2, 10000)	# Uttrykt i V_0er
R = ( sqrt( E ) - sqrt( E - 1 ) ) ** 4
T = 1 - R
plot(E, R, E, T)
legend(["R(E)", "T(E)"])
title("Plott av trans- og refleksjonskoeffisienten gitt energien til en partikkel")
xlabel("E [V_0]")
ylabel("[Aner ikke]")
savefig("Figur_05.png")
show()