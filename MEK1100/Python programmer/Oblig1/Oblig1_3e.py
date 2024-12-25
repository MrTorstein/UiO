from matplotlib.pylab import linspace, meshgrid, contour, show, savefig, title

x = linspace(-1, 1, 101)
y = linspace(-1, 1, 101)

X, Y = meshgrid(x, y)
contour(X, Y, (1 - X**2/2. - Y**2/2.))
title("Stomlinjene for (cos(x)sin(y)i - sin(x)cos(y)j), tilnermet med Taylor")
savefig("Oblig1_3e.png")
show()

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig1_3e.py

"""