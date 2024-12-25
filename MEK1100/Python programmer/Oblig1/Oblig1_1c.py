#coding=utf-8
from matplotlib.pylab import linspace, plot, show, legend, title, xlabel, ylabel, savefig, pi, tan
theta = [0, pi/5., pi/4., pi/3.]

def plotting(theta):
	t = linspace(0, 1, 100)
	x = t
	y = tan(theta)*t - tan(theta)*t**2
	plot(x, y)

plotting(theta[1])
plotting(theta[2])
plotting(theta[3])
title("(x, y)")
legend(["theta = $pi/5$", "theta = $pi/4$", "theta = $pi/3$"], loc = "best")
xlabel("x*")
ylabel("y*")
savefig("oblig1_1c.png")
show()

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig1_1c.py

"""