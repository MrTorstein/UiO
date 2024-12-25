from numpy import linspace, sqrt
from astropy.constants import c
from astropy.units import meter
from astropy.cosmology import WMAP9
from matplotlib.pyplot import figure, plot, xlabel, ylabel, savefig, show

z = linspace(0, 10, 100)
tau = 2 * c * 6.65 * 10 ** (-29) * meter ** 2 * 0.7 * meter ** (-3) / (3 * WMAP9.H(0) * 0.308) * (sqrt(0.308 * (1 + z) ** 3 + 0.692) + 1)
tau = tau.si
figure()
plot(z, tau)
xlabel(r"z")
ylabel(r"$\tau_e$")
savefig("Figure10.png")
show()