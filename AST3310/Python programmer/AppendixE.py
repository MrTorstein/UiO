""" Skrevet for Python 3 """

"""
# Likningene som trengs #

F_c = rho cp T sqrt(g delta) Hp^(-3/2) (lm/2)^2 deldel^(3/2)
"""

from sys import exit
from math import log10
from AppendixD import ModellSol
from matplotlib.pyplot import *
from scipy.optimize import fsolve, root
from numpy import zeros, ones, where, array, pi, sqrt
from AppendixC import Energy                            # Henter funksjonen for å regne ut epsilon


class ModellSol2(ModellSol):
    
    def __init__(self, L, R, M, rho, T):
        """ Definerer variabler """
        ModellSol.__init__(self, L, R, M, rho, T)
        self.alpha      = 1
        self.L_0        = 1.0 * self.L
        self.R_0        = 1.0 * self.R
        self.M_0        = 1.0 * self.M
        self.rho_0      = 1.42 * 10 ** (-7) * rho
        self.N          = 1 / (self.m_u * self.mu)
        self.expol      = False
        self.advars     = False

    def Sjekk_konvek(self, kappa, M, r, P, rho, T, L):
        g       = self.G * M / r ** 2
        self.Hp = P / (rho * g)
        lm      = self.alpha * self.Hp
        cp      = 5 / 2 * self.N * self.k
        U       = 64 * self.SB * T ** 3 * sqrt(self.Hp) / ( 3 * kappa * rho ** 2 * cp * sqrt(g) )
        
        nabad   = 2/5
        nabrad  = 3 * kappa * rho * self.Hp * L / ( 64 * self.SB * T ** 4 * pi * r ** 2 )
        ksi     = fsolve(lambda x: x ** 3 + U / lm ** 2 * x ** 2 + 4 * U ** 2 / lm ** 4 * x + U / lm ** 2 * (2 / 5 - nabrad), 1)[0]
        nabstar = nabrad - (3 * rho ** 2 * cp * sqrt(g) * self.Hp ** (- 1 / 2) * lm ** 2 * ksi ** (3) * kappa) / (64 * self.SB * T ** 3)
        nabp    = nabstar - ksi ** 2

        FC      = rho * cp * T * sqrt(g) * self.Hp ** (-3 / 2) * (lm / 2) ** 2 * ksi ** (3)

        return nabrad, nabstar, nabp, nabad, FC

    def Loopen(self):
        """ Hovedløkka til programmet, som regner ut de fem differenciallikningene ved hjelp av tidligere funksjoner. """
        d   = 10 ** 5                                           # Antall steg vi ønsker
        p   = 0.00025                                             # Parameter til bruk i variabel steglengde hvis dette tas i bruk

        # Bestemme om vi skal bruke variabel steglengde #
        a = str(input("Slå på variabel steglengde? [J/N]: "))
        if a != "J":
            dm  = - self.M_0 / d
        
        # bestemme om testverdier skal printes #
        b = str(input("Printe testverdier underveis i utregningen? [J/N]: "))
        if b == "J":
            c = int(input("Hvor mange runder mellom hver printing? [heltall]: "))
            
            # Printer tabell start for verdier underveis i løkka #
            print("|   i   |     M     |     r     |     P     |     L     |     T     |")

        # Lager array til difflikningene og legger inn initialbetingelsene #
        r       = zeros(d); r[0] = self.R_0
        P       = zeros(d); P[0] = self.get_P(self.rho_0, self.T_0)
        L       = zeros(d); L[0] = self.L_0
        T       = zeros(d); T[0] = self.T_0
        M       = zeros(d); M[0] = self.M_0
        Nabrad  = zeros(d)
        Nabstar = zeros(d)
        Nabp    = zeros(d)
        Nabad   = zeros(d)
        Fc      = zeros(d)

        self.Skaff_opacity()

        for i in range(0, d - 1):
            rho         = self.get_rho(P[i], T[i])  # Tettheten hentes
            eps         = Energy(T[i], rho).kalk()  # Energien per masse hentes
            kappa       = self.polate(T[i], rho)    # Finner opacity
            
            # Finner nablaene #
            nabrad, nabstar, nabp, nabad, FC = self.Sjekk_konvek(kappa, M[i], r[i], P[i], rho, T[i], L[i])

            # Lager en steglengde hvis variable steglengde er slått på #
            if a == "J":
                if nabad < nabp:
                    dm  = - min([abs(p * r[i] * (4 * pi * r[i] ** 2 * rho)), \
                          abs(p * P[i] * (- 4 * pi * r[i] ** 4) / (self.G * M[i])), \
                          abs(p * L[i] / eps), \
                          abs(p * T[i] * (4 * pi * r[i] ** 2 * self.Hp * rho) / (- T[i] * nabstar))])
                else:
                    dm  = - min([abs(p * r[i] * (4 * pi * r[i] ** 2 * rho)), \
                          abs(p * P[i] * (- 4 * pi * r[i] ** 4) / (self.G * M[i])), \
                          abs(p * L[i] / eps), \
                          abs(p * T[i] * (256 * pi ** 2 * self.SB * r[i] ** 4 * T[i] ** 3) / (- 3 * kappa * L[i]))])

            # Regner ut differensiallikningene #
            Nabrad[i]   = nabrad
            Nabstar[i]  = nabstar
            Nabp[i]     = nabp
            Nabad[i]    = nabad
            Fc[i]       = FC

            M[i + 1] = M[i] + dm
            r[i + 1] = r[i] + 1 / (4 * pi * r[i] ** 2 * rho) * dm
            P[i + 1] = P[i] - (self.G * M[i]) / (4 * pi * r[i] ** 4) * dm
            L[i + 1] = L[i] + eps * dm

            if nabad < nabp:
                T[i + 1] = T[i] - T[i] * nabstar / (4 * pi * r[i] ** 2 * self.Hp * rho) * dm
        
            else:
                T[i + 1] = T[i] - (3 * kappa * L[i]) / (256 * pi ** 2 * self.SB * r[i] ** 4 * T[i] ** 3) * dm

            # Printer testverdier hvis dette er slått på #
            if b == "J":
                if i // c > (i - 1) // c:
                    print("|%7i|%11.4e|%11.4e|%11.4e|%11.4e|%11.4e|"%(i + 1, M[i + 1], r[i + 1], P[i + 1], L[i + 1], T[i + 1]))
        
        # Regner ut den siste av nablaverdiene #
        Nabrad[-1], Nabstar[-1], Nabp[-1], Nabad[-1], Fc[-1] = self.Sjekk_konvek(kappa, M[-1], r[-1], P[-1], rho, T[-1], L[-1])

        # Definerer listene som klassevariabler for plotting #
        self.Nabrad     = Nabrad
        self.Nabstar    = Nabstar
        self.Nabp       = Nabp
        self.Nabad      = Nabad
        self.FC         = Fc
        self.M_         = M
        self.r_         = r
        self.P_         = P
        self.L_         = L
        self.T_         = T

        # printer siste verdi fra listene #
        print("r = %.5e"%r[-1])
        print("P = %.5e"%P[-1])
        print("L = %.5e"%L[-1])
        print("T = %.5e"%T[-1])

    def Plotter(self):
        """ Plotter verdiene som vi fikk ved å løse differenciallikningene """

        # Avgjør om figurer skal lagres og skaffer første figurnummer #
        lagre = str(input("Ønsker du å lagre figurene? [J/N]: "))
        if lagre == "J":
            fignr = int(input("Hva er første figurnummeret du skal lagre? Figurene blir lagret på formatet <Figur##.png>. [Heltall]: "))

        # Skalerer radius #
        x   = self.r_ / self.R
        y1  = self.Nabrad
        y2  = self.Nabstar
        y3  = self.Nabad
        
        # Plotter tre nablaene mot radius i samme figur #
        figure()
        title("Nabla against radius")
        plot(x, y1, x, y2, x, y3)
        xlabel("$R/R_{sun}$")
        ylabel("$nabla$")
        legend(["$nab_{stable}$", "$nab^{*}$", "$nab_{ad}$"])
        if lagre == "J":
            savefig("Figur%02i.png"%fignr)
        
        # Skaffer tettheten #
        rho = zeros(len(self.P_))
        for i in range(len(self.P_)):
            rho[i] = self.get_rho(self.P_[i], self.T_[i])

        # Skalerer verdier #
        x   = self.r_ / self.R
        y1  = self.M_ / self.M
        y2  = rho / self.rho
        y3  = self.L_ / self.L
        y4  = self.T_ / 10 ** 6
        y5  = self.P_ / self.P_[0]
        
        # Plotter alle listene mot masse i hver sin figur #
        figure()
        title("Mass against radius")
        plot(x, y1)
        xlabel("$R/R_{sun}$")
        ylabel("$M/M_{sun}$")
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 1))
        
        figure()
        title("Density against radius")
        yscale("log")
        plot(x, y2)
        xlabel("$R/R_{sun}$")
        ylabel("$rho/rho_{sun}$")
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 2))
        
        figure()
        title("Luminosity against radius")
        plot(x, y3)
        xlabel("$R/R_{sun}$")
        ylabel("$L/L_{sun}$")
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 3))
        
        figure()
        title("Temperatur against radius")
        plot(x, y4)
        xlabel("$R/R_{sun}$")
        ylabel("$T * 10^{-6}$")
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 4))
        
        figure()
        title("Pressure against radius")
        yscale("log")
        plot(x, y5)
        xlabel("$R/R_{sun}$")
        ylabel("$P/P_{sun}$")
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 5))

        # Viser figurene #
        show()

    def Cross_section(self):
        R_values    = self.r_ / self.R_0
        L_values    = self.L_ / self.L_0
        F_C_list    = self.FC
        n           = len(R_values)
        R0          = self.R_0
        show_every  = 50
        core_limit  = 0.995

        # -------------------------------------------------------------------------------------------------------
        # Assumptions:
        # -------------------------------------------------------------------------------------------------------
        # * R_values is an array of radii [unit: R_sun]
        # * L_values is an array of luminosities (at each r) [unit: L_sun]
        # * F_C_list is an array of convective flux ratios (at each r) [unit: relative value between 0 and 1]
        # * n is the number of elements in both these arrays
        # * R0 is the initial radius
        # * show_every is a variable that tells you to show every ...th step (show_every = 50 worked well for me)
        # * core_limit = 0.995 (when L drops below this value, we define that we are in the core)
        # -------------------------------------------------------------------------------------------------------

        # Cross-section of star
        # (See https://stackoverflow.com/questions/9215658/plot-a-circle-with-pyplot and
        # https://stackoverflow.com/questions/27629656/label-for-circles-in-matplotlib for details)

        figure()
        fig = gcf() # get current figure
        ax = gca()  # get current axis
        rmax = 1.2 * R0
        ax.set_xlim(-rmax, rmax)
        ax.set_ylim(-rmax, rmax)
        ax.set_aspect('equal')	# make the plot circular
        j = show_every
        for k in range(0, n - 1):
        	j += 1
        	if j >= show_every:	# don't show every step - it slows things down
        		if(L_values[k] > core_limit):	# outside core
        			if(F_C_list[k] > 0.0):		# convection
        				circR = Circle((0, 0), R_values[k], color ='red', fill=False)
        				ax.add_artist(circR)
        			else:				# radiation
        				circY = Circle((0, 0), R_values[k], color ='yellow', fill=False)
        				ax.add_artist(circY)
        		else:				# inside core
        			if(F_C_list[k] > 0.0):		# convection
        				circB = Circle((0, 0), R_values[k], color ='blue', fill = False)
        				ax.add_artist(circB)
        			else:				# radiation
        				circC = Circle((0, 0), R_values[k], color ='cyan', fill = False)
        				ax.add_artist(circC)
        		j = 0
        circR = Circle((2 * rmax, 2 * rmax), 0.1 * rmax, color = 'red', fill=True)		# These are for the legend (drawn outside the main plot)
        circY = Circle((2 * rmax, 2 * rmax), 0.1 * rmax, color = 'yellow', fill=True)
        circC = Circle((2 * rmax, 2 * rmax), 0.1 * rmax, color = 'cyan', fill=True)
        circB = Circle((2 * rmax, 2 * rmax), 0.1 * rmax, color = 'blue', fill=True)
        ax.legend([circR, circY, circC, circB], ['Convection outside core', 'Radiation outside core', 'Radiation inside core', 'Convection inside core']) # only add one (the last) circle of each colour to legend
        legend(loc = 2)
        xlabel('R')
        ylabel('R')
        title('Cross-section of star')

        # Show all plots
        show()

if __name__ == "__main__":
    L   = 3.828 * 10 ** 26
    R   = 6.957 * 10 ** 8
    M   = 1.9885 * 10 ** 30
    rho = 1.408 * 10 ** 3
    T   = 5770

    A = ModellSol2(L, R, M, rho, T)
    A.Skaff_opacity()
    A.Loopen()
    A.Plotter()
    A.Cross_section()