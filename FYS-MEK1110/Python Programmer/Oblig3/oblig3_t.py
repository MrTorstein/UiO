#coding=utf-8
from oblig3_s import txvF
from matplotlib.pylab import figure, plot, legend, xlabel, ylabel, title, savefig, show

figure()
t, x, v, F = txvF(m = 0.1, k = 100., v0 = 0, tslutt = 2., u = 0.1, n = 50000)
plot(t, F)
legend(["F(t)"])
xlabel("t")
ylabel("F")
title("Kraften paa kloss med masse 0.1kg")
savefig("oblig3_t1.png")

figure()
t_, x_, v_, F_ = txvF(m = 1.0, k = 100., v0 = 0, tslutt = 2., u = 0.1, n = 50000)
plot(t_, F_)
legend(["F(t)"], loc = "best")
xlabel("t")
ylabel("F")
title("Kraften paa kloss med masse 1.0kg")
savefig("oblig3_t2.png")

show()

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>oblig3_t.py

"""