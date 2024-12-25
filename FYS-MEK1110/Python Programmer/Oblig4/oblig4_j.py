#coding=utf-8
from oblig4_h import xt
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, title, savefig
x, t = xt(10, - 5)
plot(t, x)
legend(["x(t)"], loc = "best")
xlabel("t")
ylabel("x")
title("posisjonen til en partikkel med utgangsfart $10m/s$ i posisjon $x = -5$")
savefig("oblig4_j.png")
show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig4_j.py

"""