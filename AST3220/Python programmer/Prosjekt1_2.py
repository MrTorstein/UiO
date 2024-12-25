"""
Equations:

X^2 = Sum_{i = 1}^N ( ( dL( z_i, p ) - dL_i )^2 / sig_i^2 )
"""

from sys import exit
from numpy import array, linspace, sqrt, zeros, sum as Sum, where, ones, min as Min
from numpy.ma import masked_where
from matplotlib.pyplot import figure, plot, xlabel, ylabel, show, scatter, errorbar, savefig, legend
from Prosjekt1_1 import LumDisCalculator

"""
Constants:
X2 = Minimization quantity
"""

class Modelfinder():
    """
    Class calculating the best fit model of a LCDM (Lambda C D M) universe given measurements of redshift, luminosity distance and error, and lists of Omega_m0 and Omega_rc to try

    Takes the values:
    - List_m0   = List of Omega_m0 values to try out
    - List_L0   = Corresponding list of Omega_Lambda0 values to try out
    - File      = Datafile containing table of measured values for z, dL and error in this order. Autoset to data.txt
    
    Functions:
    - Minimize(self)
    - Find_Area(self)
    """

    def __init__(self, List_m0, List_L0, File = "data.txt"):
        """ Function to define variables as class variables """

        self.c_H0_to_Gpc    = 4.4                       # [Gpc * H_0 / c]
        
        self.List_m0        = array(List_m0)            # List of Omega_m0 values
        self.List_L0        = array(List_L0)            # List of Omega_L0 values

        self.lm0            = len(List_m0)              # Number of Omega_m0 values given
        self.lL0            = len(List_L0)              # Number of Omega_L0 values given
        if self.lm0 != self.lL0:                        # Make sure given equal numbers of Omegas
            print("Give equally many Omega_m0 and Omega_L0 values")
            exit()
        
        self.X2     = zeros([self.lm0, self.lL0])       # Xi quadrat matrix

        with open(File, "r") as infile:                 # Extracting measured data from File
            for i in range(5):                          # Skipping first 5 lines
                infile.readline()
        
            z               = []                        # List of redshift data
            dL              = []                        # List of luminosity distance data
            err             = []                        # List of errors

            for line in infile:                     # Reading data into lists
                z.append(float(line.split()[0]))
                dL.append(float(line.split()[1]))
                err.append(float(line.split()[2]))
        
        # Turning lists to class value numpy.array #
        self.z   = array(z)
        self.dL  = array(dL)
        self.err = array(err)
    
    def Minimize(self):
        """
        Function minimizing X^2 in Baye's theorem with a Gaussian distriburion
        """

        # Rename class variables to write less #
        List_m0 = self.List_m0
        List_L0 = self.List_L0
        z       = self.z
        dL      = self.dL
        err     = self.err

        # Calculate the X^2 for all the different luminosity distances for each Omega values #
        for i in range(self.lm0):
            for j in range(self.lL0):
                A = LumDisCalculator(List_m0[i], List_L0[j], False)

                self.X2[i, j] = Sum( ( A.Run(z) * self.c_H0_to_Gpc - dL ) ** 2 / err )

        self.X2_min = Min(self.X2)              # Find the minimum value in X^2
        i, j = where(self.X2 == self.X2_min)    # Find the coordinates of min_X^2
        self.Best_Omg_m0 = List_m0[i]           # Find corresponding Omega_m0 value
        self.Best_Omg_L0 = List_L0[j]           # Find corresponding Omega_L0 value

        print("Minimum value of X^2 = ", self.X2_min)
        print("Best fit Omega_m0    = ", self.Best_Omg_m0)
        print("Best fit Omega_L0    = ", self.Best_Omg_L0)

    def Find_area(self):
        """
        Finds the area around the minimum value of X2 where there is a 95% chance the true value is
        """

        # Rename class variables to write less #
        List_m0 = self.List_m0
        List_L0 = self.List_L0
        X2      = self.X2
        X2_min  = self.X2_min

        x       = []
        y       = []

        # Plot the area around the minimum and the minimum itself #
        figure()
        mask = masked_where(X2 < 6.17 + X2_min, X2).mask
        for i in range(len(X2)):
            x.append(List_m0[i] * mask[i])
            y.append(List_L0 * mask[i])
        
        return x, y

        show()

if __name__ == "__main__":

    List_m0 = linspace(0, 1, 10)
    List_L0 = linspace(0, 1, 10)

    A = Modelfinder(List_m0, List_L0)

    A.Minimize()
    
    """
    x, y = A.Find_area()
    scatter(x, y, c = "b")
    scatter(A.Best_Omg_m0, A.Best_Omg_L0, c = "r")
    xlabel("$\Omega_{m0}$")
    ylabel("$\Omega_{\Lambda 0}$")
    savefig("Figure_13.png")
    """

    figure()
    scatter(A.z, A.dL)
    errorbar(A.z, A.dL, yerr = A.err, linestyle = "None")

    B = LumDisCalculator(Omg_m0 = A.Best_Omg_m0, Omg_L0 = A.Best_Omg_L0, Printing = False)
    dL = B.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    C = LumDisCalculator(Omg_m0 = 0.2, Omg_L0 = 0.2, Printing = False)
    dL = C.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    D = LumDisCalculator(Omg_m0 = 0, Omg_L0 = 0.3, Printing = False)
    dL = D.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    E = LumDisCalculator(Omg_m0 = 0, Omg_L0 = 0.4, Printing = False)
    dL = E.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    F = LumDisCalculator(Omg_m0 = 0.3, Omg_L0 = 0.7, Printing = False)
    dL = F.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    G = LumDisCalculator(Omg_m0 = 0, Omg_L0 = 1, Printing = False)
    dL = G.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    H = LumDisCalculator(Omg_m0 = 1, Omg_L0 = 0, Printing = False)
    dL = H.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    xlabel("$z$")
    ylabel("$d_L [Gpc]$")
    legend(["$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$ (Best fit)"%(B.Omg_m0, B.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$"%(C.Omg_m0, C.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$"%(D.Omg_m0, D.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$"%(E.Omg_m0, E.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$"%(F.Omg_m0, F.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$ (dS)"%(G.Omg_m0, G.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$ (EdS)"%(H.Omg_m0, H.Omg_L0), \
            "Measured data from data.txt file" ])

    savefig("Figure_14.png")
    show()
