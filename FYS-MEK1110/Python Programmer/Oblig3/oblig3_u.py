#coding=utf-8
from oblig3_s import txvF
from matplotlib.pylab import plot, legend, xlabel, ylabel, title, savefig, show

t, x, v, F = txvF(m = 0.1, k = 10., v0 = 0, tslutt = 2., u = 0.1, n = 50000)
plot(t, x)
legend(["x(t)"], loc = "best")
xlabel("t")
ylabel("x")
title("Bevegelsen til kloss med fjerkonstant paa 10N og masse 0.1")
savefig("oblig3_u.png")

show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig3_u.py

"""