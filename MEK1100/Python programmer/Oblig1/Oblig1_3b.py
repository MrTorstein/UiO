#coding=utf-8
from matplotlib.pylab import linspace, meshgrid, savefig, sin, cos, quiver, show, title
x = linspace(-10, 10, 31)
y = linspace(-10, 10, 31)
X, Y = meshgrid(x, y)

vX = cos(X)*sin(Y)
vY = -sin(X)*cos(Y)

quiver(X, Y, vX, vY)
title("Stomvektorer for (cos(x)sin(y)i - sin(x)cos(y)j)")
savefig("oblig1_3b.png")
show()

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig1_3b.py

"""