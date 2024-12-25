#coding=utf-8
from matplotlib.pylab import linspace, meshgrid, contour, log, show, title, savefig
x = linspace(-10, 10, 500)
y = log(abs(x))
X, Y = meshgrid(x, y)
contour(X, Y, Y - log(abs(x)))
title("Stromlinjene for hastighetsfeltet v")
savefig("oblig1_2b.png")
show()
"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig1_2b.py

"""