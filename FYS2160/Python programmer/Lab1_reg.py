from numpy import polyfit, polyval, linspace, sqrt, log as ln
from matplotlib.pyplot import plot, legend, title, xlabel, ylabel, show, figure

Sample = [ [6, 735.5], [7, 863.8], [8, 1005], [9, 1144], [10, 1283], [11, 1425], [12, 1567], [13, 1709], [14, 1849], [15, 1991] ]

p = polyfit(Sample[:, 0], Sample[:, 1], 2)

x = linspace(0, 10, 1e3)
y = polyval(p, x)

#c = sqrt( (f + 2) )

figure()
plot(x, y, [1, 2], [1, 2])
title("Regresjon av resonans vs. temperatur")
legend(["Res(T)", "C(T)"], ["R", "B"])
xlabel("n [freedom]")
ylabel("f [1/s]")
show()