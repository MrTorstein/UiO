from matplotlib.pylab import plot, show, savefig, linspace, array, pi, exp, cos, sin
t = linspace(0, 4*pi, 10001)
r = lambda t: exp(-t)*array((cos(t), sin(t)))
plot(r(t)[0, :], r(t)[1, :])
savefig("Oblig1_3a.png")
show()
"""
C:\Users\Torstein\Documents\UiO\Mat1110\Python programmer>Oblig1_3a.py

"""