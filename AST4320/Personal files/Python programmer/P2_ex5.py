from numpy import linspace, e, array
from matplotlib.pyplot import figure, plot, legend, xlabel, ylabel, savefig, show

M = linspace(0.1, 2, 10 ** 5)

n_1 = 1 / M ** 2
n_2 = 1 / M ** (2/3) * e ** (-M ** (4/3))

figure()
plot(M, n_1, M, n_2)
legend([r"$M^{-2}$", r"$M^{-\frac{2}{3}} \exp(-M^{\frac{4}{3}})$"])
savefig("Figure04.png")

show()