from numpy import linspace, zeros, pi, log10, sqrt
from matplotlib.pyplot import figure, plot, xscale, yscale, ylim, xlabel, ylabel, legend, savefig, show
from astropy.constants import G
from astropy.cosmology import WMAP9 as cosmo
H = cosmo.H

a = 10 ** linspace(-3, 0, 10000)
da = a[-1] - a[0] / len(a)
adot_a = H(0) / a ** (3 / 2)

Del = zeros(len(a)); Del[0] = 10 ** (-3)
d_Del = zeros(len(a)); d_Del[0] = 3 * Del[0] * H(0) ** 2 / (2 * adot_a[0] ** 2 * (a[0] ** (-1) * Del[0] + 2 * a[0] ** (-2)))

def solve_diff(Omg_m):
    for i in range(0, len(a) - 1):
        dd_Del = - 2 * 1 / a[i] * d_Del + Del[i] * 3 / 2 * H(0) ** 2 * Omg_m / (a[i] ** 5 * adot_a[i] ** 2)
        d_Del[i + 1] = d_Del[i] + da
        Del[i + 1] = Del[i] + d_Del[i + 1] * da
    
    return d_Del, Del

figure()

d_Del_1, Del_1 = solve_diff(1)
plot(log10(a), log10(Del))

adot_a = H(0) * sqrt(0.3 / a ** (3) + 0.7)
d_Del_2, Del_2 = solve_diff(0.3)
plot(log10(a), log10(Del))

adot_a = H(0) * sqrt(0.8 / a ** (3) + 0.2)
d_Del_3, Del_3 = solve_diff(0.8)
plot(log10(a), log10(Del))

xlabel(r"$\log(a)$")
ylabel(r"$\log(\delta)$")
legend([r"$ \Omega_m = 1, \Omega_{\Lambda} = 0$", r"$\Omega_m = 0.3, \Omega_{\Lambda} = 0.7$", r"$\Omega_m = 0.8, \Omega_{\Lambda} = 0.2$"])
savefig("Figure01.png")

show()

figure()
plot(1 / a - 1, a / Del_1 * d_Del_1)
plot(1 / a - 1, a / Del_2 * d_Del_2)
plot(1 / a - 1, a / Del_3 * d_Del_3)

yscale("log")
xlabel(r"z")
ylabel(r"$f = d\log(\delta)/d\log(a)$")
legend([r"$ \Omega_m = 1, \Omega_{\Lambda} = 0$", r"$\Omega_m = 0.3, \Omega_{\Lambda} = 0.7$", r"$\Omega_m = 0.8, \Omega_{\Lambda} = 0.2$"])
savefig("Figure02.png")

show()