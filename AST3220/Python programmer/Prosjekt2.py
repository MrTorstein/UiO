from numpy import linspace, zeros, exp, log, sqrt, argmax
from numpy.ma import masked_inside
from scipy.integrate import ode
from matplotlib.pyplot import rcParams, figure, plot, xlabel, ylabel, ylim, xscale, yscale, legend, savefig, show

class WIMPModel():

    def __init__(self, m_, cross, steps):
        """
        Definerer konstanter og variabler

        """
        S   = 1 / 2
        self.m_     = m_
        self.cross  = cross
        self.g      = 2 * S + 1
        self.g_     = 7 / 8 * 2 * self.g
        self.dx     = 1000 - 1 / steps
        
        self.x      = linspace(1, 1000, steps)
        self.W      = zeros(steps)

        self.y1     = zeros(1)
        self.y1_eq  = zeros(1)
        self.y2     = zeros(1)
        self.y2_eq  = zeros(1)
        self.y3     = zeros(1)
        self.y3_eq  = zeros(1)

        # Spesifiserer større font størrelse for figurene #
        self.fs     = 14
        rcParams.update({"font.size": self.fs})

    def _Loopen(self):
        """
        Integrasjonsloopen
        """

        self.W_eq   = log(9.35 * 10 ** 9 * self.g / 2 * sqrt(100 / self.g_) * self.m_ / 1000 * self.cross * 10 ** (10) * self.x ** (3 / 2) * exp(-self.x))

        self.W[0] = self.W_eq[0]
        
        f = lambda x, W, W_eq: x ** (-2) * (exp(2 * W_eq - W) - exp(W))

        temp = ode(f).set_integrator("vode", method = "bdf")
        temp.set_initial_value(self.W[0], self.x[0])
        
        for i in range(len(self.x) - 1):
            temp.set_f_params(self.W_eq[i])
            self.W[i + 1] = temp.integrate(self.x[i + 1])

    def _Plotter(self):
        """
        Plotter resultatene
        """
        figure()
        plot(self.x, self.y1, color = "red")
        plot(self.x, self.y1_eq, color = "green")
        plot(self.x, self.y2, color = "blue")
        plot(self.x, self.y2_eq, color = "turquoise")
        plot(self.x, self.y3, color = "orange")
        plot(self.x, self.y3_eq, color = "purple")

        ylim((0.5, 1.5 * 10 ** 12))
        xscale("log")
        yscale("log")
        xlabel("x")
        ylabel("y")
        legend([r"$y(x, \langle \sigma v \rangle = 10^{-9})$", r"$y_{eq}(x, \langle \sigma v \rangle = 10^{-9})$", \
                r"$y(x, \langle \sigma v \rangle = 10^{-10})$", r"$y_{eq}(x, \langle \sigma v \rangle = 10^{-10})$", \
                r"$y(x, \langle \sigma v \rangle = 10^{-11})$", r"$y_{eq}(x, \langle \sigma v \rangle = 10^{-11})$"])
        
        savefig("Figure_20.png")
        show()
    
    def Run(self):
        """
        Kjører størsteparten av modellen for tre forskjellige cross
        """
        self._Loopen()
        self.y1     = exp(self.W)
        self.y1_eq  = exp(self.W_eq)

        self.cross *= 10 ** (-1)
        self._Loopen()
        self.y2     = exp(self.W)
        self.y2_eq  = exp(self.W_eq)

        self.cross *= 10 ** (-1)
        self._Loopen()
        self.y3     = exp(self.W)
        self.y3_eq  = exp(self.W_eq)

        self._Plotter()

    def CalOmg(self, cross):
        """
        Regner ut mengde mørk materie idag
        """
        x_f = self.x[argmax(self.y1 < 0.1 * self.y1[0])]
        return 1.69 * x_f / 20 * sqrt(100 / self.g_) * 10 ** (-10) / cross

if __name__ == "__main__":
    A = WIMPModel(1000, 10 ** (-9), 10 ** 4)
    A.Run()

    Omega1      = A.CalOmg(10 ** (-9))
    Omega2      = A.CalOmg(10 ** (-10))
    Omega3      = A.CalOmg(10 ** (-11))
    print("Omega h^2(<ov> = 10^(-9))  = %.3f" %Omega1)
    print("Omega h^2(<ov> = 10^(-10)) = %.3f" %Omega2)
    print("Omega h^2(<ov> = 10^(-11)) = %.3f" %Omega3)

    crosses     = linspace(10 ** (-14), 10 ** (-7), 10 ** 8)
    Omg         = A.CalOmg(crosses)
    interval    = crosses[crosses * masked_inside(Omg, 0.07, 0.17).mask > 0]
    print("<ov> = %.3e +- %.3e" %(interval[0] + (interval[-1] - interval[0]) / 2, (interval[-1] - interval[0]) / 2))
