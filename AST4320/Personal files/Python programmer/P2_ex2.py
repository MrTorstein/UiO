from numpy import linspace, sin, pi
from matplotlib.pyplot import figure, plot, xlabel, ylabel, legend, savefig, show

R = pi

k = linspace(-20, 20, 1000)
W = 2 * sin(R * k) / k

figure()
plot(k, W)
xlabel("k")
ylabel(r"$\tilde{W}(k)$")
legend([r"$R = \pi$"])
savefig("Figure03.png")
show()