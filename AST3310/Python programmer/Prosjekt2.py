""" Skrevet for Python 3 """

"""
# Likningene som trengs #

dr / dm = 1 / (4 pi r^2 rho)
dP / dm = - (G m) / (4 pi r^4)
dL / dm = eps
P       = P_G + P_rad
dT / dm = - nabla T G m / (4 pi P r^4)
P_G     = rho k T / (mu m_u)
P_rad   = a T^4 / 3
Nabla   = dln(T) / dln(P)
"""

from sys import exit
from math import log10
from numpy import zeros, ones, where, array, pi, sum as Sum, sqrt, cbrt, roots, isreal, array
from Prosjekt1 import Energy                          # Henter funksjonen for å regne ut epsilon
from matplotlib.pylab import plot, show, figure, title, savefig, xlabel, ylabel, yscale, legend, ylim, xlim, rcParams
from cross_section import cross_section

class ModellSol():
    """
    Klasse for å modelere energitransporten i den sentrale delen av en sol liknende stjerne.
    Tar variablene:
    - L                     = Raten av solas luminositeten modelen skal bygges på 
    - R                     = Raten av solas radiusen modelen skal bygges på
    - M                     = Raten av solas massen modelen skal bygges på
    - rho                   = Raten av solas gjennomsnittstettheten modelen skal bygges på
    - T                     = Temperaturen model skal bygges på
    - d                     = Maksimalt antall steg for å begge modellen
    - sanity                = variabel for å slå på og av 

    Inneholder funksjonene:
    - _Skaff_opacity(self)
    - _polate(self, T, rho, sanity = False)
    -- T            = temperatur
    -- rho          = massetetthet
    -- sanity       = parameter gitt i funksjonen sanitytest for å teste interpoleringsmetoden
    - _get_p(self, rho, T)
    -- rho          = massetetthet
    -- T            = temperatur
    - _get_rho(self, P, T)
    -- P            = totalt trykk
    -- T            = temperatur
    - _dRdm(self, rho, R)
    -- rho          = massetetthet
    -- R            = radius
    - _dPdm(self, M, R)
    -- M            = masse
    -- R            = radius
    - _dLdm(self, rho, T)
    -- rho          = massetetthet
    -- T            = temperatur
    - _dTdm(self, M, T, P, R)
    -- M            = masse
    -- T            = temperatur
    -- P            = totalt trykk
    -- R            = radius
    - _nabla_stable(self, rho, L, M, T, R, kappa)
    -- rho          = massetetthet
    -- L            = luminositet
    -- M            = masse
    -- T            = temperatur
    -- R            = radius
    -- kappa        = opasitet
    - _nabla_ad(self)
    - _nabla_star(self, rho, M, T, R, kappa)
    -- rho          = massetetthet
    -- M            = masse
    -- T            = temperatur
    -- R            = radius
    -- kappa        = opasitet
    - _nabla_p(self)
    - _Loopen(self, p = 0.01, a = "N", b = "N")
    -- p            = parameter for å bestemme største mulige massesteget ved variabel steglengde
    -- a            = parameter for å slå på og av variabel steglengde 
    -- b            = parameter for å slå på og av printing av parametere underveis
    - _Cross_section_plotter(self, show_every = 1, nr = 10, lagre = "N")
    -- show_every   = parameter for å bestemme mellomrommet mellom hvert lag i modellen som skal plottes
    -- nr           = tall brukt for å bestemme navn på lagret figur
    -- lagre        = parameter for å bestemme om figurer skal lagres
    -- fs           = fontsize of plot
    - _Sanitytest(self)
    - Run(self)
    - Plotter(self)


    Kjøres ved å kalle på Run. Deretter kan du hente ut variablene:
    - self.L_               = d elementer lang liste med luminositetsverdier
    - self.R_               = d elementer lang liste med radiusverdier
    - self.M_               = d elementer lang liste med massesverdier
    - self.P_               = d elementer lang liste med trykkverdier
    - self.T_               = d elementer lang liste med temperaturverdier
    - self.nabla_stable_    = d elementer lang liste med temperatur gradient verdier
    - self.nabla_star_      = d elementer lang liste med temperatur gradient verdier
    - self.nabla_ad_        = d elementer lang liste med temperatur gradient verdier
    
    Resultater plottes ved Plotter()
    """
    
    def __init__(self, L, R, M, rho, T, P = False, d = 10 ** 4, sanity = None):
        """ Definerer variabler """

        self.L          = L                                 # Ratio av solas luminositet
        self.R          = R                                 # Ratio av solas radius
        self.M          = M                                 # Ratio av solas masse
        self.rho        = rho                               # Ratio av solas gjennomsnittstetthet
        self.P          = P                                 # Ratio av solas fotosfæriske gass trykk
        self.L_S        = 3.846 * 10 ** 26                  # Luminositeten til sola [W]
        self.R_S        = 6.96 * 10 ** 8                    # Radien til sola [m]
        self.M_S        = 1.989 * 10 ** 30                  # Massen til sola [kg]
        self.rho_S      = 1.408 * 10 ** 3                   # Gjennomsnitstettheten til sola [kg/m^3]
        self.P_S        = 10 ** 4                           # Fotosfærisk gass trykk [N/m^2]
        self.G          = 6.674 * 10 ** (-11)               # Newtons gravitasjonskonstant [m^3kg^-1s^-2]
        self.m_u        = 1.660539040 * 10 ** (-27)         # Atomisk masseenhet [kg]
        self.SB         = 5.670367 * 10 ** (-8)             # Stefan-Boltzmanns konstant [W/m^2K^4]
        self.c          = 299792458                         # Lyshastigheten [m/s]
        self.N_A        = 6.022140857 * 10 ** 23            # Avogadros konstant [1/mol]
        self.k          = 1.38064852 * 10 ** (-23)          # Boltzmanns konstant [J/K]
        self.pi         = 3.141592654                       # Pi
        self.a          = 4 * self.SB / self.c              # Strålingstetthet konstant [J/m^3K^4]
        self.L_0        = self.L * self.L_S                 # Luminositeten til modellen vår [W]
        self.R_0        = self.R * self.R_S                 # Radius til modellen vår [m]
        self.M_0        = self.M * self.M_S                 # Utgangsmassen til hele vår modell [kg]
        self.rho_0      = self.rho * self.rho_S             # Utgangstettheten til vår modell [kg/m^3]
        self.T_0        = T                                 # Utgangstemperaturen til modellen vår [K]
        self.P_0        = self.P * self.P_S                 # Utgangstrykket til modellen vår [N/m^2]
        self.alpha      = 1                                 # Parameter mellom 1/2 og 2 for å bestemme mikselengden
        self.X          = 0.7                               # Hvor stor del av modellen vår sin masse som kommer av hydrogen
        self.Y          = 0.29                              # Hvor stor del av modellen vår sin masse som kommer av helium 4
        self.Y_3        = 10 ** (-10)                       # Hvor stor del av modellen vår sin masse som kommer av helium 3
        self.Z_Li       = 10 ** (-7)                        # Hvor stor del av modellen vår sin masse som kommer av litium
        self.Z_Be       = 10 ** (-7)                        # Hvor stor del av modellen vår sin masse som kommer av berylium
        self.Z_N        = 10 ** (-11)                       # Hvor stor del av modellen vår sin masse som kommer av nitrogen
        self.mu         = 1 / (2 * self.X + 3 / 4 * self.Y + self.Y_3 + 4 / 7 * self.Z_Li + 5 / 7 * self.Z_Be + 8 / 14 * self.Z_N) # Gjennomsnittlig molekylvekt
        self.c_p        = 5 * self.k / (2 * self.mu * self.m_u)                                                                    # Spesifik varmekapasitet

        self.d          = d                                 # Antall steg i modellen vår

        # Lager array til difflikningene og legger inn initialbetingelsene for modelen #
        self.R_ = zeros(d + 1); self.R_[0] = self.R_0
        self.L_ = zeros(d + 1); self.L_[0] = self.L_0
        self.T_ = zeros(d + 1); self.T_[0] = self.T_0
        self.M_ = zeros(d + 1); self.M_[0] = self.M_0
        self.P_ = zeros(d + 1)
        if self.P != False:
            self.P_[0] = self.P_0
        else:
            self.P_[0] = self._get_P(self.rho_0, self.T_0)

        
        # Arrays for å lagre nabla verdier #
        self.nabla_stable_  = zeros(d + 1)
        self.nabla_ad_      = zeros(d + 1)
        self.nabla_star_    = zeros(d + 1)
        
        # Parametere for å advare om extrapolering #
        self.expol      = False
        self.advars     = False

        # Spørr om sanitytest skal slås av eller på hvis dette ikke er prebestemt #
        while sanity != "J" and sanity != "N":
            sanity = input("Sanity test slått på? [J/N]: ")
        
        self.sanity = sanity

        # Spesifiserer større font størrelse for figurene #
        self.fs     = 22
        #rcParams.update({"font.size": self.fs})
        
    def _Skaff_opacity(self):
        """ 
        Leser opacity.txt og henter ut verdiene for opasiteten som en numpy.array. 
        Henter også de tilsvarende verdiene for LogT og LogR i arrays 
        """

        # Finner dimmensjonene til matrisa først #
        with open("opacity.txt", "r") as innfil:
            counter1 = 0
            counter2 = 0

            innfil.readline()
            innfil.readline()
            for line in innfil:
                counter1 += 1
                counter2 = len(line.split()[1:])
            dim = [counter1, counter2]
        
        # Henter ut første rad, første rekke og selve kappa matrisa #
        with open("opacity.txt", "r") as innfil:
            LogR = innfil.readline().split()[1:]
            LogR = ones(len(LogR)) * array( list( map(float, LogR) ) )
            Kappa = zeros(dim)
            LogT = zeros(dim[0])
            innfil.readline()
            i = 0
            for line in innfil:
                LogT[i] = line.split()[0]
                Kappa[i] = line.split()[1:]
                i += 1
        
        # Definerer som klassevariabler for senere bruk #
        self.LogR   = LogR
        self.LogT   = LogT
        self.Kappa  = Kappa

    def _polate(self, T, rho, sanity = False):
        """
        Interpolerer, eller extrapolerer for å finne opasitetsverdien ved T og rho

        Trenger at self._Skaff_opacity() er kjørt før bruk 

        Tar verdiene:
        - T         = temperatur
        - rho       = massetetthet
        - sanity    = parameter gitt i funksjonen sanitytest for å teste interpoleringsmetoden
        
        Returnerer opasitetsverdien ved temperaturen T og massetettheten rho
        """
        
        # Sjekker om funksjonen blir kalt for testing eller for bruk #
        if sanity == False:
            x       = log10(T)
            y       = log10(rho / (1000 * (T * 10 ** (-6)) ** 3)) # skaffer log10(R) verdien i cgs
        
        else:
            x = T
            y = rho
        
        # Henter listene som poleres over #
        LogR    = self.LogR
        LogT    = self.LogT
        Kappa   = self.Kappa

        # Definerer lister over tall som brukes i poleringen #
        b = zeros(2)
        c = zeros(2)

        # Finner to verdier av temperaturne å bruke #
        if x > LogT[-1]:                            # Ekstrapolerer hvis T er for stor
            self.expol  = True
            tall1       = -2
            b[0]        = LogT[-2]
            b[1]        = LogT[-1]
        elif x < LogT[0]:                           # Ekstrapolerer hvis T er for liten
            self.expol  = True
            tall1       = 0
            b[0]        = LogT[1]
            b[1]        = LogT[0]
        else:                                       # Interpolerer
            a = LogT[LogT > x][0]
            tall1   = where(a == LogT)[0][0] - 1
            b[0]    = a
            b[1]    = LogT[tall1]

        # Finner to verdier av R å bruke #
        if y > LogR[-1]:                            # Ekstrapolerer
            self.expol  = True
            tall2       = -2
            c[0]        = LogR[-2]
            c[1]        = LogR[-1]
        elif y < LogR[0]:                           # Ekstrapolerer
            self.expol  = True
            tall2       = 0
            c[0]        = LogR[1]
            c[1]        = LogR[0]
        else:                                       # Interpolerer
            a = LogR[LogR > y][0]
            tall2 = where(a == LogR)[0][0] - 1
            c[0] = LogR[tall2]
            c[1] = a

        # Bruker interpoleringsmetode hentet fra wikipedia til å finne en tilnærmet verdi #
        f = ( Kappa[tall1 + 1, tall2] * (c[1] - y) * (b[1] - x) \
            + Kappa[tall1 + 1, tall2 + 1] * (y - c[0]) * (b[1] - x) \
            + Kappa[tall1, tall2] * (c[1] - y) * (x - b[0]) \
            + Kappa[tall1, tall2 + 1] * (y - c[0]) * (x - b[0]) ) \
              / ( (c[1] - c[0]) * (b[1] - b[0]) )

        # Printer en advarsel hvis deet blir ekstrapolert #
        if self.expol == True and self.advars == False:
            print("ADVARSEL!!! Du har nå extrapolert utenfor de tilgjengelige verdiene av opasiteten og dette er en ganske dårlig tilnærming")
            self.advars = True

        return 10 ** f * 0.1                        # Gjør poleringsverdien om til SI enheter

    def _get_P(self, rho, T):
        """ 
        Finner det totale trykket, inni modellen vår, ved hjelp av ideel gass lov og strålingstrykk
        
        Tar verdiene:
        - rho   = massetetthet
        - T     = temperatur

        Returnerer det totale trykket
        """

        return rho * self.k * T / (self.mu * self.m_u) + self.a * T ** 4 / 3

    def _get_rho(self, P, T):
        """
        Finner tettheten ved hjelp av ideel gass lov og strålingstrykk
        
        Tar verdiene:
        - P = totalt trykk
        - T = temperatur
        
        Returnerer massetettheten
        """

        return (P - self.a * T ** 4 / 3) * self.mu * self.m_u / ( self.k * T )

    def _dRdm(self, rho, R):
        """
        Finner endringen i radius, til modellen, per endring i masse

        Tar verdiene:
        - rho   = massetetthet
        - R     = radius

        Returnerer endringen i radius per endring i masse
        """

        return 1 / (4 *  self.pi * R ** 2 * rho)

    def _dPdm(self, M, R):
        """
        Finner endringen i totalt trykk, til modellen, per endring i masse
        
        Tar verdiene:
        - M = masse
        - R = radius

        Returnerer endringen i det totalet trykket per endring i masse
        """
        
        return - (self.G * M) / (4 * self.pi * R ** 4)

    def _dLdm(self, rho, T):
        """
        Finner endringen i luminositet, til modellen, per endring i masse
        
        Tar verdiene:
        - rho   = massetetthet
        - T     = temperatur

        Returnerer endringen i luminositet per endring i masse
        """

        return Energy(T, rho, "N").Kalk()

    def _dTdm(self, M, T, P, R):
        """
        Finner endringen i temperatur, til modellen, per endring i masse
        
        Krever at self._nabla_star(rho, M, T, R, kappa) er kjørt først for å gi korrekt verdi

        Tar verdiene:
        - M = masse
        - T = temperatur
        - P = totalt trykk
        - R = radius

        Returnerer endringen i temperatur per endring i masse
        """

        return self.nabla_star * T / P * self._dPdm(M, R)

    def _nabla_stable(self, rho, L, M, T, R, kappa):
        """
        Regner ut både nabla_stable og H_p og setter dem som klassevariabler
        
        Tar verdiene:
        - rho   = massetettheten
        - L     = luminositeten
        - M     = massen
        - T     = temperaturen
        - R     = radius
        - kappa = opasiteten

        Definerer trykkskaleringshøyden og temperaturgradienten ved kun strålingstransport som klassevariabler
        """

        self.H_p            = - self._get_P(rho, T) * self._dRdm(rho, R) / self._dPdm(M, R)
        self.nabla_stable   = 3 * kappa * rho * self.H_p * L / (64 * self.pi * R ** 2 * self.SB * T ** 4)

    def _nabla_ad(self):
        """ Regner ut nabla_ad og setter den som klassevariabel """
        self.nabla_ad = 2 / 5

    def _nabla_star(self, rho, M, T, R, kappa):
        """
        Regner ut g, U, l_m, xi og nabla_star, og setter dem som klassevariabler

        Krever at self._nabla_stable(rho, L, M, T, R, kappa) er kjørt først for å gi riktig svar
        
        Tar verdiene:
        - rho   = massetetthet
        - M     = masse
        - T     = temperatur
        - R     = radius
        - kappa = opasitet

        Definerer gravitasjonsakselrasjonen, konstanten U, blandingslengden, konstanten xi og verdien for temperaturgradienten ved et lag i modellen 
        som klassevariabler
        """

        # Regner ut variabler som trengs #
        g                   = self.G * M / R ** 2 
        U                   = 64 * self.SB * T ** 3 / (3 * kappa * rho ** 2 * self.c_p) * sqrt(self.H_p / g)
        l_m                 = self.alpha * self.H_p

        # Sjekker for konveksjon og regner ut nabla_star #
        if self.nabla_stable < self.nabla_ad:
            self.nabla_star = self.nabla_stable
        
        else:
            xi = roots([1, U / (l_m ** 2), 4 * U ** 2 / (l_m ** 4), U / (l_m ** 2) * (self.nabla_ad - self.nabla_stable)])
            xi = xi[isreal(xi)].real[0]
        
            """
            # Analytisk løsning av xi. Upraktisk å bruke, men ser pen ut da =) #
            xi                  = - 11 * cbrt(2) * U ** 2 \
                                  / (3 * cbrt(- 27 * nabla_ad * l_m ** 10 * U \
                                              + sqrt(((- 27 * nabla_ad * l_m ** 10 * U \
                                                                  + 27 * l_m ** 10 * nabla_stable * U \
                                                                  + 34 * l_m ** 6 * U ** 3)) ** 2 \
                                                                  + (5324 * l_m ** 12 * U ** 6)) \
                                              + 27 * l_m ** 10 * nabla_stable * U \
                                              + 34 * l_m ** 6 * U ** 3)) \
                                  + cbrt(- 27 * nabla_ad * l_m ** 10 * U \
                                         + sqrt(((- 27 * nabla_ad * l_m ** 10 * U \
                                                             + 27 * l_m ** 10 * nabla_stable * U \
                                                             + 34 * l_m ** 6 * U ** 3)) ** 2 \
                                                             + (5324 * l_m ** 12 * U ** 6)) \
                                         + 27 * l_m ** 10 * nabla_stable * U \
                                         + 34 * l_m ** 6 * U ** 3) \
                                  / (3 * cbrt(2) * l_m ** 4) \
                                  - U / (3 * l_m ** 2)
            """
            
            # Definerer variabler som klassevariabler #
            self.xi             = xi
            self.nabla_star     = self.nabla_stable - l_m ** 2 / U * xi ** 3

        # Definerer variabler som klassevariabler #
        self.g              = g
        self.U              = U
        self.l_m            = l_m

    def _nabla_p(self):
        """
        Regner ut nabla_p og setter den som klassevariabel
        
        Krever at self._nabla_star(rho, M, T, R, kappa) er kjørt først for å gi riktig svar
        """
        xi = roots([1, self.U / (self.l_m ** 2), 4 * self.U ** 2 / (self.l_m ** 4), self.U / (self.l_m ** 2) * (self.nabla_ad - self.nabla_stable)])
        xi = xi[isreal(xi)].real[0]
        self.nabla_p = 4 * self.U / self.l_m ** 2 * xi + self.nabla_ad

    def _Loopen(self, p = 0.01, a = "N", b = "N"):
        """
        Hovedløkka til programmet, som regner ut de fire differensiallikningene.
        
        Tar verdiene:
        - p = parameter for å bestemme største mulige massesteget ved variabel steglengde
        - a = parameter for å slå på og av variabel steglengde
        - b = parameter for å slå på og av printing av parametere underveis

        Definerer listene L_, R_, M_, P_, T_, nabla_stable_, nabla_star_, nabla_ad_ som klassevariabler
        """

        # Finner steglengden hvis den er konstant #
        dm  = - 0.99 * self.M_0 / self.d

        # Setter opp tabell hvis dette skal printes #
        if b == "J":
            c = int(input("Hvor mange runder mellom hver printing? [Heltall]: "))
            print("|  i  |     R     |     M     |     P     |     L     |     T     |Nabla_stable| Nabla_star |    rho    |")
        
        i = 0
        
        while self.M_[i] + dm > 0 and self.R_[i] > 1 and self.P_[i] > 1 and self.L_[i] > 1 and self.T_[i] > 1:
            rho     = self._get_rho(self.P_[i], self.T_[i])  # Finner tettheten
            kappa   = self._polate(self.T_[i], rho)          # Henter opasiteten
            
            # Regner ut nabla verdier #
            self._nabla_stable(rho, self.L_[i], self.M_[i], self.T_[i], self.R_[i], kappa)  # Finner nabla_stable
            self._nabla_ad() 
            self._nabla_star(rho, self.M_[i], self.T_[i], self.R_[i], kappa)
            
            # Lagrer nabla verdier #
            self.nabla_stable_[i]   = self.nabla_stable
            self.nabla_ad_[i]       = self.nabla_ad
            self.nabla_star_[i]     = self.nabla_star
            
            # Finner stigningen til differensiallikningene #
            f = [ \
                 self._dRdm(rho, self.R_[i]), \
                 self._dPdm(self.M_[i], self.R_[i]), \
                 self._dLdm(rho, self.T_[i]), \
                 self._dTdm(self.M_[i], self.T_[i], self.P_[i], self.R_[i])]
            
            # Lager en steglengde hvis variable steglengde er slått på #
            if a == "J":
                dm  = - min([ \
                      abs(p * self.R_[i] / f[0]), \
                      abs(p * self.P_[i] / f[1]), \
                      abs(p * self.L_[i] / f[2]), \
                      abs(p * self.T_[i] / f[3])])
            
            # Regner ut differensiallikningene #
            self.M_[i + 1] = self.M_[i] + dm
            self.R_[i + 1] = self.R_[i] + f[0] * dm
            self.P_[i + 1] = self.P_[i] + f[1] * dm
            self.L_[i + 1] = self.L_[i] + f[2] * dm
            self.T_[i + 1] = self.T_[i] + f[3] * dm
            
            # Printer testverdier hvis dette er slått på #
            if b == "J":
                if i // c > (i - 1) // c:
                    print("|%5i|%11.4e|%11.4e|%11.4e|%11.4e|%11.4e|%12.4e|%12.4e|%11.4e|"%(i, self.R_[i], self.M_[i], self.P_[i], self.L_[i], \
                                                                                           self.T_[i], self.nabla_stable, self.nabla_star, rho))
            
            # Øker indeksen med en #
            i += 1
        
        # Sletter plasser i lister som ikke er brukt #
        self.M_ = self.M_[self.M_ > 0]
        self.R_ = self.R_[self.R_ > 0]
        self.P_ = self.P_[self.P_ > 0]
        self.L_ = self.L_[self.L_ > 0]
        self.T_ = self.T_[self.T_ > 0]
        self.nabla_stable_  = self.nabla_stable_[0:len(self.M_)]
        self.nabla_ad_      = self.nabla_ad_[0:len(self.M_)]
        self.nabla_star_    = self.nabla_star_[0:len(self.M_)]

        # Regner ut de aller siste nabla verdiene #
        M       = self.M_[-1]
        P       = self.P_[-1]
        T       = self.T_[-1]
        R       = self.R_[-1]
        L       = self.L_[-1]
        rho     = self._get_rho(P, T)
        kappa   = self._polate(T, rho)

        self._nabla_stable(rho, L, M, T, R, kappa)
        self._nabla_ad()
        self._nabla_star(rho, M, T, R, kappa)

        self.nabla_stable_[-1]  = self.nabla_stable
        self.nabla_star_[-1]    = self.nabla_star
        self.nabla_ad_[-1]      = self.nabla_ad

    def _Cross_section_plotter(self, show_every = 1, nr = 10, lagre = "N", fs = 10):
        """
        Plotter et utsnitt av stjerna ved å kalle på funksjonen fra cross_section.py
        
        Tar verdiene:
        - show_every    = parameter for å bestemme mellomrommet mellom hvert lag i modellen som skal plottes
        - nr            = tall brukt for å bestemme navn på lagret figur
        - lagre         = parameter for å bestemme om figurer skal lagres
        - fs            = fontsize of plot
        """

        # Regner ut verdier som trengs, og som ikke normalt lagres #
        rho_    = self._get_rho(self.P_, self.T_)
        H_p_     = - self.P_ * self._dRdm(rho_, self.R_) / self._dPdm(self.M_, self.R_)
        kappa_  = zeros(len(self.M_))
        for i in range(len(kappa_)):
            kappa_[i]   = self._polate(self.T_[i], rho_[i])
        F_c = 16 * self.SB * self.T_ ** 4 * (self.nabla_stable_ - self.nabla_star_) / (3 * kappa_ * rho_ * H_p_)

        # Kjører cross_section.py #
        cross_section(self.R_, self.L_, F_c, show_every, nr, lagre, fs)

    def _Sanitytest(self):
        """ Gjennomfører tester for klassen """

        ## Test en ##
        print("## Test en ##")
        # Definerer tilfeldige testvariabler #
        T       = 10 ** 6
        rho     = 10 ** 3
        P       = 10

        # Printer testverdier for å sjekke om P(rho(P, T), T) = P og rho(P(rho, T), T) = rho #
        print("rho' = rho(P(rho', T), T)")
        print("%d = %d"%(rho, self._get_rho(self._get_P(rho, T), T)))
        print("P' = P(rho(P', T), T)")
        print("%d = %d"%(P, self._get_P(self._get_rho(P, T), T)))

        print("-----------------------------------------------------------------")


        ## Test to ##
        print("## Test to ##")
        # Sjekker om interpoleringen funker #
        logTer  = [3.750, 3.755, 3.755, 3.755, 3.755, 3.770, 3.780, 3.795, 3.770, 3.775, 3.780, 3.795, 3.800]
        logRer  = [-6.00, -5.95, -5.80, -5.70, -5.55, -5.95, -5.95, -5.95, -5.80, -5.75, -5.70, -5.55, -5.50]
        self._Skaff_opacity()
        print("|%6s|%6s|%8s|"%("log(T)", "log(R)", "kappa"))
        for i in range(len(logTer)):
            print("|%6.3f|%6.3f|%8.2e|"%(logTer[i], logRer[i], self._polate(logTer[i], logRer[i], True)))
            
        print("-----------------------------------------------------------------")


        ## Test tre ##
        print("## Test tre ##")
        # Definerer testvariabler #
        T                   = 0.9 * 10 ** 6                 # [T]
        rho                 = 55.9                          # [kg m^-3]
        R                   = 0.84 * self.R_S               # [m]
        M                   = 0.99 * self.M_S               # [kg]
        kappa               = 3.98                          # [m^2 kg^-1]
        alpha               = 1

        # Regner ut variabler som trengs #
        L                   = self.L_0
        self._nabla_stable(rho, L, M, T, R, kappa)
        self._nabla_ad()
        self._nabla_star(rho, M, T, R, kappa)
        v                   = sqrt(self.g * self.l_m ** 2 / (4 * self.H_p)) * self.xi
        F_c                 = (self.nabla_stable - self.nabla_star) / self.nabla_stable
        F_r                 = self.nabla_star / self.nabla_stable
        self._nabla_p()

        # Plotter testverdier for sammenlikning #
        print("Nabla_stable = %.2f"%self.nabla_stable)
        print("Nabla_ad     = %.1f"%self.nabla_ad)
        print("H_p          = %.2e"%self.H_p)
        print("U            = %.2e"%self.U)
        print("xi           = %.3e"%self.xi)
        print("Nabla_star   = %.6f"%self.nabla_star)
        print("v            = %.2f"%v)
        print("F_c/F_cr     = %.2f"%F_c)
        print("F_r/F_cr     = %.2f"%F_r)
        print("Nabla_p      = %.6f"%self.nabla_p)

        print("-----------------------------------------------------------------")


        ## Test fire ##
        print("## Test fire ##")
        # Definerer initial parametere #
        self.L_[0]   = self.L_S
        self.R_[0]   = self.R_S
        self.M_[0]   = self.M_S
        self.T_[0]   = 5770
        self.P_[0]   = self._get_P(self.rho_S * 1.42 * 10 ** (-7), self.T_[0])

        # Finner Nabla verdiene #
        self._Loopen(a = "J",b = "J")

        # Plotter Nabla verdier #
        print("Plotter Nabla verdier")
        figure(figsize = (9.5, 7))
        plot(self.R_ / self.R_S, self.nabla_stable_)
        plot(self.R_ / self.R_S, self.nabla_star_)
        plot(self.R_ / self.R_S, self.nabla_ad_)
        legend([r"$\nabla_{stable}$", r"$\nabla^{*}$", r"$\nabla_{ad}$",], fontsize = self.fs)
        xlabel(r"$R/R_{\odot}$", fontsize = self.fs)
        ylabel(r"$\nabla$", fontsize = self.fs)
        yscale("log")
        ylim((0.1, 10 ** 3))
        savefig("Figur00.png")
        show()
        
        print("-----------------------------------------------------------------")


        ## Test fem ##
        print("## Test fem ##")
        print("Plotter snitt av stjerna")
        self._Cross_section_plotter(5, lagre = "J", fs = self.fs)

        print("-----------------------------------------------------------------")


        # Avslutter når tester er ferdig #
        exit()

    def Run(self):
        """ Funksjon for å kjøre modellen """
        if self.sanity == "J":
            self._Sanitytest()
        
        p   = 0.01                          # Parameter til bruk i variabel steglengde hvis dette tas i bruk

        # Bestemme om vi skal bruke variabel steglengde #
        a = str(input("Slå på variabel steglengde? [J/N]: "))
        
        # Bestemme om testverdier skal printes og hvor ofte #
        b = str(input("Printe testverdier underveis i utregningen? [J/N]: "))
        
        # Skaffer opasitetsmatrisa #
        self._Skaff_opacity()

        # Løser differensiallikningene #
        self._Loopen(p, a, b)

        # printer siste verdi fra listene #
        print("R  = %.5e"%self.R_[-1])
        print("P  = %.5e"%self.P_[-1])
        print("L  = %.5e"%self.L_[-1])
        print("T  = %.5e"%self.T_[-1])

    def Plotter(self):
        """
        Plotter verdiene som vi fikk ved å løse differensiallikningene
        Trenger at Run() er kjørt først
        """

        # Avgjør om figurer skal lagres og skaffer første figurnummer #
        lagre = str(input("Ønsker du å lagre figurene? [J/N]: "))
        if lagre == "J":
            fignr = int(input("Hva er første figurnummeret du skal lagre? Figurene blir lagret på formatet <Figur##.png>. [Heltall]: "))
        else:
            fignr = 80

        # Skaffer tettheten #
        rho = self._get_rho(self.P_, self.T_)

        # Skalerer verdiene #
        x   = self.R_ / self.R_S
        y1  = self.M_ / self.M_S
        y2  = self.T_ / 10 ** 6
        y3  = self.L_ / self.L_S
        y4  = rho
        y5  = self.P_
        
        # Plotter alle listene mot radius i hver sin figur #
        figure()
        plot(x, y1)
        xlabel(r"$R/R_{\odot}$" )
        ylabel(r"$M/M_{\odot}$")
        legend([r"$M(R)$"])
        if lagre == "J":
            savefig("Figur%02i.pdf"%fignr)
        
        figure(figsize = (9.5, 7))
        plot(x, y2)
        xlabel(r"$R/R_{\odot}$", fontsize = self.fs)
        ylabel(r"$T \cdot 10^{-6} [K]$", fontsize = self.fs)
        legend([r"$T(R)$"], fontsize = self.fs)
        ylim((-0.3, 18.3))
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 1))
        
        figure(figsize = (8, 7))
        plot(x, y3)
        xlabel(r"$R/R_{\odot}$", fontsize = self.fs)
        ylabel(r"$L/L_{\odot}$", fontsize = self.fs)
        legend([r"$L(R)$"], fontsize = self.fs)
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 2))
        
        figure(figsize = (10, 7))
        plot(x, y4)
        xlabel(r"$R/R_{\odot}$", fontsize = self.fs)
        ylabel(r"$\rho$ [kg/m^3]", fontsize = self.fs)
        legend([r"$\rho(R)$"], fontsize = self.fs)
        yscale("log")
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 3))
        
        figure(figsize = (10, 7))
        plot(x, y5)
        xlabel(r"$R/R_{\odot}$", fontsize = self.fs)
        ylabel(r"$P [N/m^2]$", fontsize = self.fs)
        legend([r"$P(R)$"], fontsize = self.fs)
        yscale("log")
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 4))

        show() # Viser figurene

        # Finner strålingsfluks- og konveksjonsfluksraten av total fluks #
        F_c = (self.nabla_stable_ - self.nabla_star_) / self.nabla_stable_
        F_r = self.nabla_star_ / self.nabla_stable_

        # Plotter Fluks mot radius #
        figure(figsize = (8, 7))
        plot(x, F_c, x, F_r)
        xlabel(r"$R/R_{\odot}$", fontsize = self.fs)
        ylabel(r"Flux ratio", fontsize = self.fs)
        legend([r"$F_{c}/F_{tot}$", r"$F_{r}/F_{tot}$"], fontsize = self.fs)
        savefig("Figur%02i.png"%(fignr + 5))

        show()

        # Finner energiutbyttet fra hver av kjedene #
        eps0 = []
        eps1 = []
        eps2 = []
        eps3 = []
        eps4 = []
        for i in range(len(self.T_)):
            B = Energy(self.T_[i], rho[i], "N")
            B.Kalk()
            eps0.append(B.eps[0])
            eps1.append(B.eps[1])
            eps2.append(B.eps[2])
            eps3.append(B.eps[3])
            eps4.append(B.eps[4])
        eps0 = array(eps0)
        eps1 = array(eps1)
        eps2 = array(eps2)
        eps3 = array(eps3)
        eps4 = array(eps4)
        
        # Plotter energiutbytte mot radius #
        figure(figsize = (8.5, 7.5))
        plot(x * 10 ** 4, eps0 / eps4)
        plot(x * 10 ** 4, eps1 / eps4)
        plot(x * 10 ** 4, eps2 / eps4)
        plot(x * 10 ** 4, eps3 / eps4)
        plot(x * 10 ** 4, eps4 / max(eps4))
        xlabel(r"$R/R_{\odot} \cdot 10^{-4}$", fontsize = self.fs)
        ylabel(r"$\epsilon [Js^{-1}kg^{-1}]$", fontsize = self.fs)
        legend([r"$\epsilon_{PPI}/\epsilon(R)$", \
                r"$\epsilon_{PPII}/\epsilon(R)$", \
                r"$\epsilon_{PPIII}/\epsilon(R)$", \
                r"$\epsilon_{CNO}/\epsilon(R)$", \
                r"$\epsilon(R)/\epsilon_{max}$"], fontsize = self.fs)
        xlim((0, 3))
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 6))
        
        show()

        # Plotter Nabla verdier #
        figure(figsize = (9.5, 7))
        plot(self.R_ / self.R_S, self.nabla_stable_)
        plot(self.R_ / self.R_S, self.nabla_star_)
        plot(self.R_ / self.R_S, self.nabla_ad_)
        legend([r"$\nabla_{stable}$", r"$\nabla^{*}$", r"$\nabla_{ad}$"], fontsize = self.fs)
        xlabel(r"$R/R_{\odot}$", fontsize = self.fs)
        ylabel(r"$\nabla$", fontsize = self.fs)
        yscale("log")
        ylim((0.1, 10 ** 3))
        if lagre == "J":
            savefig("Figur%02i.png"%(fignr + 7))
        show()

        # Plotter snittet av stjerna #
        if input("Plott utsnitt av stjerna? [J/N]: ") == "J":
            self._Cross_section_plotter(5, fignr + 8, lagre, self.fs)


if __name__ == "__main__":
    L   = 1.0
    R   = 1.0 * 1.02
    M   = 1.0
    rho = 1.42 * 10 ** (-7) * 20
    T   = 5770 * 0.2
    P   = False

    A = ModellSol(L, R, M, rho, T, P, d = 7 * 10 ** 6)
    A.Run()
    A.Plotter()
