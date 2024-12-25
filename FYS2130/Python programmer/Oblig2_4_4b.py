from Oblig2_4_4 import RK4, plotte
from matplotlib.pylab import plot, legend, title, xlabel, ylabel, savefig, show

"""
Variabler:
"""
b 	= 10	#kg/s
k 	= 10	#N/m
m 	= 0.1	#kg
z_0	= 0.1	#m
v_0	= 0		#m/s
T	= 10	#s
N	= 1e4

t1, z1, v1 = RK4(- b / m, - k / m, 0, z_0, v_0, T, N)

b = 2

t2, z2, v2 = RK4(- b / m, - k / m, 0, z_0, v_0, T, N)

b = 0.1

t3, z3, v3 = RK4(- b / m, - k / m, 0, z_0, v_0, T, N)


plot(t1, z1, t2, z2, t3, z3)
title("Funksjon med %d antall tidssteg" %(N))
xlabel("tid [s]")
ylabel("posisjon [m]")
legend(["Overkritisk b = 10", "Kritisk b = 2", "Underkritisk b = 0.1"])
savefig("Figur_0004.png")
show()
