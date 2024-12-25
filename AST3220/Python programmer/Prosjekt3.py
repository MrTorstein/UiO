from scipy.constants import pi
from numpy import zeros, sqrt, sum, append
from matplotlib.pyplot import figure, plot, xscale, yscale, xlabel, ylabel, legend, savefig, show

class Slow_Roller():
    """

    """
    def __init__(self):
        self.n      = 10 ** 5
        self.dtau   = 2.5 * 10 ** 2 / self.n

        self.tau    = zeros(self.n); self.tau[0]    = 0
        self.h      = zeros(self.n); self.h[0]      = 1
        self.psi    = zeros(self.n); self.psi[0]    = 3.1
        self.psi_v  = zeros(self.n)
        self.lna_ai = zeros(self.n); self.lna_ai[0] = self.h[0] / 2 * self.dtau
    
    def _Loopen(self):
        self.psi_v[0] = - 1 / (4 * pi * self.h[0] * self.psi[0])
        for i in range(self.n - 1):
            dv_dpsi             = 3 * self.psi[i] / (4 * pi * self.psi[0] ** 2)
            psi_a               = - 3 * self.h[i] * self.psi_v[i] - dv_dpsi
            self.psi_v[i + 1]   = self.psi_v[i] + psi_a * self.dtau

            self.tau[i + 1]     = self.tau[i] + self.dtau
            self.psi[i + 1]     = self.psi[i] + self.psi_v[i + 1] * self.dtau
            self.h[i + 1]       = sqrt(8 * pi / 3 * (1 / 2 * self.psi_v[i + 1] ** 2 + 3 * self.psi[i + 1] ** 2 / (8 * pi * self.psi[0] ** 2)))

            self.lna_ai[i + 1]  = self.lna_ai[i] + self.h[i + 1] * self.dtau
        
        self.lna_ai[-1] = self.lna_ai[-2] + self.h[-1] / 2 * self.dtau
    
    def _Plotter(self):
        figure()
        plot(self.tau, self.psi)
        plot(self.tau, self.psi[0] - 1 / (4 * pi * self.psi[0]) * self.tau)
        xlabel(r"$\tau$")
        ylabel(r"$\psi$")
        legend(["Numerical solution", "Analytical solution"])
        savefig("Figure_30.png")

        figure()
        plot(self.tau, self.lna_ai)
        xlabel(r"$\tau$")
        ylabel(r"$\ln(a / a_i)$")
        savefig("Figure_31.png")

        v = 3 * self.psi ** 2 / (8 * pi * self.psi[0] ** 2)

        figure()
        plot(self.tau, 1 / 2 * self.psi_v ** 2 - v / (1 / 2 * self.psi_v ** 2 + v))
        xlabel(r"$\tau$")
        ylabel(r"$p_{\phi} / \rho_{\phi} c^2$")
        savefig("Figure_32.png")

        show()

    def Run(self):
        self._Loopen()
        self._Plotter()

if __name__ == "__main__":
    A = Slow_Roller()
    A.Run()