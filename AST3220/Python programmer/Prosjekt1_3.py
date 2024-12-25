"""
Equations:

FI: H^2 / H_0^2 = ( sqrt( Omg_m0 ( z + 1 )^3 + Omg_rc ) + sqrt(Omg_rc) )^2 + ( 1 - Omg_m0 - Omg_L0 ) ( z + 1 )^2
dL = c( z + 1 ) / ( H_0 sqrt( |Omega_k0 ) ) S_k( sqrt( |Omega_k0 ) * int_0^z H_0 / H dz )

X^2 = Sum_{i = 1}^N ( ( dL( z_i, p ) - dL_i )^2 / sig_i^2 )
"""

from sys import exit
from numpy import array, linspace, sqrt, sin, sinh, zeros, sum as Sum, where, ones, min as Min
from numpy.ma import masked_where
from scipy.integrate import simps
from Prosjekt1_1 import LumDisCalculator
from Prosjekt1_2 import Modelfinder

"""
Constants:
- H_0   = Hubbles constant
- h     = Dimentionless Hubbel constant
- c     = Speed of light
- X2    = Minimization quantity
"""

class LumDisCalculator2(LumDisCalculator):
    """
    Class to calculate luminosity distance at a redshift z for the DGP model.

    It takes the values:
    - Omg_m0    = Omega_m0, mass density parameter
    - Omg_rc    = Omega_Lambdrc, length scale density parameter
    - Printing  = Boolean determining if function prints or doesn't print the successful values for Omg_m0 and Omg_L0

    Functions contained:
    - _FI(self, z)
    -- z        = Redshift
    - _Test_Omega(self)
    - _Calculate_Distance(self, z, steps = 10 ** 4)
    -- z        = Redshift
    -- steps    = Number of integration steps, automatically set to 10^4
    - Run(self, zlist, steps = "False")
    -- zlist    = a z value, or a numpy.array of z values
    -- steps    = steps sent to _Calculate_Distance function. If false, auto value used
    """

    def __init__(self, Omg_m0 = None, Omg_rc = None, Printing = True):
        """ Defining constants and variables as class variables """
        
        LumDisCalculator.__init__(self, Omg_m0 = Omg_m0, Omg_L0 = 0, Printing = Printing)

        # Finding values through input #
        if Omg_rc == None:   # Obtaining Omega_Lambda0
            Omg_rc  = float(input("Give a value for Omega_rc: "))

        self.Omg_rc = Omg_rc
        self.Omg_k0 = 1 - (sqrt(Omg_m0 + Omg_rc) + sqrt(Omg_rc)) ** 2
        
    def _FI(self, z):
        """
        Calculates the right hand side of the Friedmann equation

        Takes the value:
        - z     = redshift
        """

        return ( sqrt( self.Omg_m0 * (z + 1) ** (3) + self.Omg_rc ) + sqrt(self.Omg_rc) ) ** 2 + self.Omg_k0 * (z + 1) ** (2)

    def _Test_Omega(self):
        """ Testing if Right side of Friedmann equation is positive """

        afrac = linspace(10 ** (-5), 1, 10 ** 5)    # list of a / a_0 fraction values

        RS = self._FI(1 / afrac - 1)                # Right side of FI
        RS = RS[RS < 0]                             # Remove all values above or equal to 0

        timer = 0                                   # Attempts used to find different values

        while len(RS) > 0 and timer < 5 and self.Printing: # Ask for new value and test these new values
            print("\n### Given values made Friedmanns first equation problematic. Give other values ###")
            
            self.Omg_m0 = float(input("Give a value for Omega_m0: "))
            self.Omg_rc = float(input("Give a value for Omega_rc: "))
            self.Omg_k0 = 1 - (sqrt(self.Omg_m0 + self.Omg_rc) + sqrt(self.Omg_rc)) ** 2

            RS      = self._FI(1 / afrac - 1)
            RS      = RS[RS < 0]

            timer += 1  # Add 1 to timer

        if timer == 5:  # Used to many attempts. Exit program
            print("\n### Took to long to find resonable values. Please try again =) ###")
            exit()

        elif self.Printing:
            print("You were successful with the values \n Omega_m0 = %f \n Omega_rc = %f"%(Omg_m0, self.Omg_rc))
        

class Modelfinder2(Modelfinder):
    """
    Class calculating the best fit model of a DGP universe given measurements of redshift, luminosity distance and error, and lists of Omega_m0 and Omega_rc to try

    Takes the values:
    - List_m0   = List of Omega_m0 values to try out
    - List_L0   = Corresponding list of Omega_Lambda0 values to try out
    - File      = Datafile containing table of measured values for z, dL and error in this order. Autoset to data.txt
    
    Functions:
    - Minimize(self)
    - Find_Area(self)
    """

    def __init__(self, List_m0, List_rc, File = "data.txt"):
        """ Function to define variables as class variables """

        Modelfinder.__init__(self, List_m0 = List_m0, List_L0 = List_rc, File = File)
        
        self.List_rc = List_rc  # List of Omega_rc values
    
    def Minimize(self):
        """
        Function minimizing X^2 in Baye's theorem with a Gaussian distriburion
        """

        # Rename class variables to write less #
        List_m0 = self.List_m0
        List_rc = self.List_rc
        z       = self.z
        dL      = self.dL
        err     = self.err

        # Calculate the X^2 for all the different luminosity distances for each Omega values #
        for i in range(self.lm0):
            for j in range(self.lL0):
                A = LumDisCalculator2(List_m0[i], List_rc[j], False)

                self.X2[i, j] = Sum( ( A.Run(z) * self.c_H0_to_Gpc - dL ) ** 2 / err )

        self.X2_min = Min(self.X2)              # Find the minimum value in X^2
        i, j = where(self.X2 == self.X2_min) # Find the coordinates of min_X^2
        self.Best_Omg_m0 = List_m0[i]           # Find corresponding Omega_m0 value
        self.Best_Omg_rc = List_rc[j]           # Find corresponding Omega_L0 value

        print("Minimum value of X^2 = ", self.X2_min)
        print("Best fit Omega_m0    = ", self.Best_Omg_m0)
        print("Best fit Omega_rc    = ", self.Best_Omg_rc)

if __name__ == "__main__":    
    """ """
    from matplotlib.pyplot import figure, plot, xlabel, ylabel, show, scatter, errorbar, savefig, legend

    List_m0 = linspace(0, 1, 10)
    List_L0 = linspace(0, 1, 10)
    List_rc = linspace(0, 0.3, 10)

    A = Modelfinder2(List_m0, List_rc)
    A.Minimize()
    """
    x, y = A.Find_area()
    scatter(x, y, c = "b")
    scatter(A.Best_Omg_m0, A.Best_Omg_rc, c = "r")
    xlabel("$\Omega_{m0}$")
    ylabel("$\Omega_{rc}$")
    savefig("Figure_16.png")
    """
    B = Modelfinder(List_m0, List_L0)
    B.Minimize()

    figure()
    scatter(A.z, A.dL)
    errorbar(A.z, A.dL, yerr = A.err, linestyle = "None")

    C = LumDisCalculator2(Omg_m0 = A.Best_Omg_m0, Omg_rc = A.Best_Omg_rc, Printing = False)
    dL = C.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    D = LumDisCalculator(Omg_m0 = B.Best_Omg_m0, Omg_L0 = B.Best_Omg_L0, Printing = False)
    dL = D.Run(B.z) * B.c_H0_to_Gpc
    plot(B.z, dL)

    E = LumDisCalculator(Omg_m0 = 0, Omg_L0 = 1, Printing = False)
    dL = E.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    F = LumDisCalculator(Omg_m0 = 1, Omg_L0 = 0, Printing = False)
    dL = F.Run(A.z) * A.c_H0_to_Gpc
    plot(A.z, dL)

    xlabel("$z$")
    ylabel("$d_L [Gpc]$")
    legend(["$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$ (Best fit LCDM)"%(D.Omg_m0, D.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{rc} = %2.2f$ (Best fit DGP)"%(C.Omg_m0, C.Omg_rc), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$ (dS)"%(E.Omg_m0, E.Omg_L0), \
            "$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$ (EdS)"%(F.Omg_m0, F.Omg_L0), \
            "Measured data from data.txt file" ])

    savefig("Figure_17.png")
    show()
    
