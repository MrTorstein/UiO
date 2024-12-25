"""
Equations:

FI: H^2 / H_0^2 = Omg_m0 ( z + 1 )^3 + ( 1 - Omg_m0 - Omg_L0 ) ( z + 1 )^2 + Omg_L0
dL = c( z + 1 ) / ( H_0 sqrt( |Omega_k0 ) ) S_k( sqrt( |Omega_k0 ) * int_0^z H_0 / H dz )
"""

from sys import exit
from numpy import linspace, sqrt, zeros, sin, sinh, array
from scipy.integrate import simps

"""
Constants:
- H_0     = Hubbles constant
- h       = Dimentionless Hubbel constant
- c       = Speed of light
"""
class LumDisCalculator():
    """
    Class to calculate luminosity distance at a redshift z for the LCDM model (Lambda C D M model).

    It takes the values:
    - Omg_m0    = Omega_m0, mass density parameter
    - Omg_L0    = Omega_Lambda0, cosmological constant density parameter
    - Printing  = Boolean determining if function prints or doesn't print the values for Omg_m0 and Omg_L0

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

    def __init__(self, Omg_m0 = None, Omg_L0 = None, Printing = True):
        """ Defining constants and variables as class variables """

        # Finding values through input #
        if Omg_m0 == None:   # Obtaining Omega_m0
            Omg_m0  = float(input("Give a value for Omega_m0: "))
        if Omg_L0 == None:   # Obtaining Omega_Lambda0
            Omg_L0  = float(input("Give a value for Omega_Lambda0: "))
        
        # Defining variables as class variables #
        self.Printing = Printing
        self.Omg_m0 = Omg_m0
        self.Omg_L0 = Omg_L0
        self.Omg_k0 = 1 - Omg_m0 - Omg_L0
        self.H_0    = 100 * 0.7 * 3.154 * 10 ** (7) / (3.086 * 10 ** (19))  # Hubbles constant [yr^-1]
        self.c      = 3 * 10 ** 8                                           # Speed of light [m/s]

    def _FI(self, z):
        """
        Calculates the right hand side of the Friedmann equation

        Takes the value:
        - z     = redshift
        """
        return self.Omg_m0 * (z + 1) ** (3) + ( 1 - self.Omg_m0 - self.Omg_L0 ) * (z + 1) ** (2) + self.Omg_L0

    def _Test_Omega(self):
        """ Testing if Right side of Friedmann equation is positive """
        afrac = linspace(10 ** (-5), 1, 10 ** 5)            # list of a / a_0 fraction values

        RS = self._FI(1 / afrac - 1)                        # Right side of FI
        RS = RS[RS < 0]                                     # Remove all values above or equal to 0

        timer = 0                                           # Attempts used to find different values

        while len(RS) > 0 and timer < 5 and self.Printing:  # Ask for new value and test these new values
            print("\n### Given values made Friedmanns first equation problematic. Give other values ###")
            
            self.Omg_m0 = float(input("Give a value for Omega_m0: "))
            self.Omg_L0 = float(input("Give a value for Omega_Lambda0: "))
            self.Omg_k0 = 1 - self.Omg_m0 - self.Omg_L0

            RS      = self._FI(1 / afrac - 1)
            RS      = RS[RS < 0]

            timer += 1  # Add 1 to timer

        if timer == 5:  # Used to many attempts. Exit program
            print("\n### Took to long to find resonable values. Please try again =) ###")
            exit()

        elif self.Printing:
            print("You were successful with the values \n Omega_m0 = %f \n Omega_Lambda0 = %f"%(self.Omg_m0, self.Omg_L0))
        
    def _Calculate_Distance(self, z, steps = 10 ** 4):
        """
        Function calculating the luminosity distance given a redshift. Intended as an inclass function

        Takes the values
        - z       = Redshift value we want to find the distance to
        - steps   = Number of integration steps, automatically set to 10^4
        """

        Omg_L0  = self.Omg_L0
        Omg_m0  = self.Omg_m0
        Omg_k0  = self.Omg_k0
        H_0     = self.H_0
        dz_     = z / steps
        temp    = zeros(steps); temp[0] = dz_

        # Finding S_k function and sign in front of Omega_k0 #
        sign = 1
        if Omg_k0 < 0:
            S = sin
            sign = -1
        elif Omg_k0 <= 10 ** -15:
            Omg_k0 = 10 ** (-15)
            S = lambda r: r
        else:
            S = sinh

        # Integrate expression inside S_k function #
        Z = linspace(0, z, steps)
        f = lambda z_: 1 / sqrt(self._FI(z_))

        int_f = simps(f(Z), Z)

        # Calculating distance #
        self.dL = (z + 1) / sqrt(sign * Omg_k0) * S(sqrt(sign * Omg_k0) * int_f)

        return self.dL
    
    def Run(self, zlist, steps = "False"):
        """
        Function to run the __Calculate_Distance function
        Takes the values:
        - zlist = a z value, or a numpy.array of z values
        - steps = steps sent to __Calculate_Distance function. If false, auto value used

        ## Recommend not running for more then 1000 values of z ##
        """

        # Running Omega test #
        self._Test_Omega()

        if steps == "False":                                # If steps is not given, dont give this value to __Calc
            if type(zlist) == int or type(zlist) == float:  # Run if given a number
                return self._Calculate_Distance(zlist)
            elif type(zlist) == type(zeros(2)):             # Run if given list of numbers
                dLlist = []
                for z in zlist:
                    dLlist.append(self._Calculate_Distance(z))
                return array(dLlist)
            else:                                           # Exit if given something that is neither a number or a list
                print("The variable 'zlist' must be a numpy array, int or float")
                exit()
        
        elif type(steps) == int or type(steps) == float:    # If steps is given correctly, give it to __Calc
            if type(zlist) == int or type(zlist) == float:  # Run if given a number
                return self._Calculate_Distance(zlist, steps)
            elif type(zlist) == type(zeros(2)):             # Run if given list of numbers
                dLlist = []
                for z in zlist:
                    dLlist.append(self._Calculate_Distance(z, steps))
                return array(dLlist)
            else:                                           # Exit if given something that is neither a number or a list
                print("The variable 'zlist' must be a numpy array, int or float")
                exit()
        else:                                               # Terminate program if steps is given wrongly
            print("The variable 'steps' must be an int or a float")
            exit()


if __name__ == "__main__":
    """ Generating pictures for Omega_k0 > 0, = 0 and < 0 """
    from matplotlib.pyplot import figure, plot, xlabel, ylabel, legend, savefig, show

    z = linspace(10 ** (-4), 10, 10 ** 4)
    
    A = LumDisCalculator(Omg_m0 = 0, Omg_L0 = 0)
    res = A.Run(z)
    
    figure()
    plot(z, res)
    xlabel("$z$")
    ylabel("$dL [c / H_0]$")
    legend(["$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$"%(A.Omg_m0, A.Omg_L0)])
    savefig("Figure_10.png")

    A.Omg_m0 = 1
    A.Omg_L0 = 0
    res = A.Run(z)

    figure()
    plot(z, res)
    xlabel("$z$")
    ylabel("$dL [c / H_0]$")
    legend(["$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$"%(A.Omg_m0, A.Omg_L0)])
    savefig("Figure_11.png")
    
    A.Omg_m0 = 0.1
    A.Omg_L0 = 0.1
    res = A.Run(z)

    figure()
    plot(z, res)
    xlabel("$z$")
    ylabel("$dL [c / H_0]$")
    legend(["$\Omega_{m0} = %2.2f, \Omega_{\lambda 0} = %2.2f$"%(A.Omg_m0, A.Omg_L0)])
    savefig("Figure_12.png")

    show()

    """ Testing some spesific cases that can be calculated empiricaly """
    
    z = 1
    A = LumDisCalculator(Omg_m0 = 0, Omg_L0 = 0)
    print("dL =", A.Run(z))

    z = 1
    A = LumDisCalculator(Omg_m0 = 1, Omg_L0 = 0)
    print("dL =", A.Run(z))