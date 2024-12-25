""" Egen Kode """
from ast2000tools.space_mission import SpaceMission
from Definer_fysiske_konstanter import DFK
from numpy import array, zeros, linspace, random, ma, sum as Sum, power, sqrt
from sys import exit, stdout

class Rocket_simulator():
    """
    Klasse for å simulere en rakettmotor og beregne drivstoff som skal brukes for å nå unnslippningshastighet.

    Tar imot verdiene:
    - Seed: et tall som bestemmer hva de tilfeldige verdiene er. Så det er enklere å teste funksjoner

    Inneholder funksjonene:
    - __init__(seed)
        Definerer forskjellige konstanter og slikt
    - _Sett_opp_motor()
        Setter opp vektorer for posisjon og hastighet i 3D for partikler, og gir partiklene verdier.
    - _Veggspretting()
        Sjekker om partikler treffer en vegg i motoren og korrigerer partiklers hastighet hvis de treffer
    - _Teller()
        Teller partikler som unnslipper og forskjellige andre fysiske enheter som skal lagres under simulasjonen.
    - _Lastefelt(i)
        Funksjon som setter opp og lager ett lastefelt for en for-løkke
    - _Tester()
        Sjekker enkle fysiske konsepter som skal stemme og gir tilbakemelding hvis noe er galt
    - Kjør_motorsimulasjon()
        Selve funksjonen som gjennomfører simulasjonen av motoren og gir fra seg antall partikler som sendes ut av motoren, 
        i tillegg til den totale bevegelsesmengden til partiklene som forsvinner
    - Regn_ut()
        Regner ut antall motorbokser som trengs for å nå unnslipningshastighet gitt en tid, og hvor mye drivstoff man trenger.
    - Utskytningsimulator()
        En simulator for å skyte opp en satelitt i rommet
    """

    def __init__(self, seed):
        self.seed       = seed; random.seed(seed)                       # Setter et seed for numpy.random klassen

        self.Var        = DFK(["pi", "G", "u", "k", "M_sol", "R_sol"])  # Definering av forskjellige fysiske konstanter
        Mission         = SpaceMission(self.seed)                       # Gjør det mulig å hente konstanter fra ast2000 klassen
        self.M_p        = Mission.system.masses[0] * self.Var.M_sol     # [kg] Massen til hjemplaneten
        self.R_p        = Mission.system.radii[0] * 10 ** 3             # [m] Radien til hjemplaneten
                
        self.L          = 10 ** -6                                      # [m] Lengde og bredde på boksen
        self.T          = 10 ** 4                                       # [K] Temperatur i boksen
        self.m_H        = 2 * 1.00794 * self.Var.u                      # [kg] Masse til et hydrogenmolekyl
        self.n          = 10 ** 5                                       # Antall partikler i boksen
        self.Dt         = 10 ** -9                                      # [s] Total tid for simulasjonen
        self.N          = 10 ** 3                                       # Antall steg i simulasjonen
        self.dt         = self.Dt / self.N                              # [s] Endring i tid for et steg
        self.n_fri      = 0                                             # Antall molekyler som unnslipper boksen
        self.p          = 0                                             # Totale mengden bevegelsesmengde som boksen får

        self.M          = 1000                                          # [kg] Vekta til raketten

    def _Sett_opp_motor(self):
        sigma       = sqrt(self.Var.k * self.T / self.m_H)      # Avvik i normalfordelingen til hastigheten

        self.pos    = random.random((self.n, 3)) * self.L
        self.hast   = random.normal(scale = sigma, size = (self.n, 3))

    def _Veggspretting(self):
        a = ma.masked_outside(self.pos, 0, self.L)              # Finner de partiklene som er uttafor boksen
        self.hast += -2 * self.hast * a.mask                    # Endrer hastigheten til de partiklene som er utenfor boksen med 180 grader

        b = ma.masked_inside(self.pos[:, 0], self.L / 4, 3 * self.L / 4)      # Finner alle partikler som treffer hullet
        c = ma.masked_inside(self.pos[:, 1], self.L / 4, 3 * self.L / 4)
        e = ma.masked_outside(self.pos[:, 2], 0, self.L)
        #print(Sum(b.mask*c.mask * e.mask))
        self.hullteller = Sum(e.mask * b.mask * c.mask)
        self.bevegelsesmengde = abs(Sum(self.m_H * self.hast[:, 0] * b.mask * c.mask * e.mask))


    def _Teller(self):
        self.n_fri += self.hullteller                           # Teller opp den totale mengden partikler som forsvinner fra motoren
        self.p += self.bevegelsesmengde                         # Summerer den totale bevegelsesmengden til partiklene som forvinner

    def _Lastefelt(self, i):
        if i == 0:
            print("Simulerer motor")
            bredde = 20											# Bredden på lastingsmåleren
            stdout.write("Laster... [%s]" % (" " * (bredde - 1)))
            stdout.flush()
            stdout.write("\b" * bredde) 						# Returnerer til starten av linja, etter "["
            self.mellomrom = self.N / bredde
            self.teller = 0        

        if i // self.mellomrom > self.teller:
            stdout.write("=")
            stdout.flush()
            self.teller = i / self.mellomrom

        elif i == self.N - 1:
            stdout.write("] Done!\n")

    def _Tester(self):
        gjen_abs_hast = Sum(sqrt(Sum(power(self.hast, 2), 1))) / self.n


        gjen_kin_energi = 0.5 * self.m_H * gjen_abs_hast ** 2
        analytisk = 3 / 2 * self.Var.k * self.T
        if 0.1 <= abs(gjen_kin_energi - analytisk) / analytisk <= 10:
            print("Gjennomsnittlig kinetisk energi svarer med analytisk verdi")
        else:
            print("Gjennomsnittlig kinetisk energi svarer ikke med analytisk verdi")
            print("Nummerisk = ", gjen_kin_energi)
            print("Analytisk = ", analytisk)
            exit()


        konst = self.m_H / (2 * self.Var.k * self.T)
        analytisk = (konst ** (3 / 2) / self.Var.pi) * 4 * self.Var.pi / (2 * konst ** 2)
        if 0.1 <= abs(gjen_abs_hast - analytisk) / analytisk <= 10:
            print("Gjennomsnittlig absolutt hastighet svarer med analytisk verdi")
        else:
            print("Gjennomsnittlig absolutt hastighet svarer ikke med analytisk verdi")
            print("Nummerisk = ", gjen_abs_hast)
            print("Analytisk = ", analytisk)
            exit()

    def Kjør_motorsimulasjon(self):
        self._Sett_opp_motor()
        self._Tester()

        for i in range(self.N):
            self._Lastefelt(i)

            self.pos += self.hast * self.dt

            self._Veggspretting()
            self._Teller()
        
        self._Tester()

    def Regn_ut(self):
        v_ush       = sqrt(2 * self.Var.G * self.M_p / self.R_p)                # Regner ut unnslippningshastighet
        print("Unnslipningshastigheten til planeten er: %e [m/s]" %v_ush)

        v_boks      = self.p / (self.M)                                         # Hastighet per boks
        print("Hastigheten raketten får fra en boks over en tid på 10^(-9) sek: %e [m/s]" %v_boks)
        
        n_bokser    = v_ush * self.Dt / ( v_boks * 20 * 60 )                    # Minste antall bokser
        print("Det minste antall bokser for at raketten skal nå unnslipningshastighet i løpet av 20 min: %e" %n_bokser)

        m_drivstoff = self.m_H * n_bokser * self.n_fri * 20 * 60 / self.Dt      # Drivstoff brukt
        print("Drivstoff brukt for å nå unnslipningshastighet: %e [kg]" %m_drivstoff)

if __name__ == "__main__":
    A = Rocket_simulator(12827)
    A.Kjør_motorsimulasjon()
    A.Regn_ut()