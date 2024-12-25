#coding=utf-8
from AST2000SolarSystem import AST2000SolarSystem
from matplotlib.pylab import array

seed = 12827
MittSolSystem = AST2000SolarSystem(seed)
#MittSolSystem.part2A_4(7, friend_seed = 15622)


"""
FRAME 2:

Initial posisjon til satelitt 1 = 0						km
Initial posisjon til sattelitt 2 = Avstand mellom
Avstand mellom skipene er L = 457.49 + 90.23 = 547.72	km
Hastighet til satelittene = 0.891						c
T_A = 0		millisek = 0								km
T_B = 7.97	millisek = 2391.000							km
T_C = 8.87	millisek = 2661.000							km
T_D = 16.69	millisek = 5007.000							km


Pos_sat_1(t) = Init_pos_sat1 + Hast_sat*t
Pos_sat_2(t) = Init_pos_sat2 + Hast_sat*t
Definert for t >= T_A:
Pos_beam_1(t) = Init_pos_sat1 + Hast_sat*T_A + c*(t - T_A)
Definert for t >= T_B:
Pos_beam_2(t) = Init_pos_sat2 + Hast_sat*T_B - c*(t - T_B)


Antar t >= T_B:
Pos_sat_1(t) = Pos_beam_2(t)
Init_pos_sat1 + Hast_sat*t = Init_pos_sat2 + Hast_sat*T_B - c*(t - T_B)
Hast_sat*t + c*t = Init_pos_sat2 - Init_pos_sat1 + Hast_sat*T_B + c*T_B
(Hast_sat + c)*t = Init_pos_sat2 + Hast_sat*T_B + c*T_B
	Init_pos_sat2 + Hast_sat*T_B + c*T_B
t = ------------------------------------
		      (Hast_sat + c)
t = 0.0089354856337 sek = 2680.645						km

Antar t >= T_A:
Pos_sat_2(t) = Pos_beam_1(t)
Init_pos_sat2 + Hast_sat*t = Init_pos_sat1 + Hast_sat*T_A + c*(t - T_A)
Hast_sat*t - c*t = Init_pos_sat1 + Hast_sat*T_A - Init_pos_sat2 - c*T_A
(Hast_sat - c)*t = - Init_pos_sat2
	 - Init_pos_sat2 
t = ----------------
	 (Hast_sat - c)
t = 0.0167498470948 sek = 5024.954						km


____________________________________________
| Event | Tid[ms] | Tid[km] | Posisjon[km] |
|-------+---------+---------+--------------|
|   A   |  00.00  | 0000.00 |   0000.000   |
|   B   |  07.97  | 2391.00 |   2677.381   |
|   C   |  08.94  | 2680.65 |   2388.455   |
|   D   |  16.75  | 5024.95 |   5024.954   |
|   E   |  25.69  | 7705.60 |   7413.409   |


Bruker eventene D og E:

ds_DE^2 = dt_DE^2 - dx_DE^2 = (7705.60 - 5024.95)^2 - (7413.409 - 5024.954)^2 = 1481167.13548 # burde vært rundt 4320000
ds_DE'^2 = dt_DE'^2 - dx_DE'^2 = dt^2
dt = sqrt(1481167.13548) = 1217.03210125 så hvis L' = 1217.03210125:
dt = L'
T_E er da lik 2L'

Dette kunne gjettes ut fra at lyset måtte gått fra skip 2 til skip 1 og tilbake igjen, \\
som vil si at tiden det tar for lyset å bevege seg avstanden mellom skipene to ganger, \\
og med tid oppgitt i km så er dette det samme som den avstanden.


FRAME 1:
____________________________________________
| Event | Tid[ms] | Tid[km] | Posisjon[km] |
|-------+---------+---------+--------------|
|   A   |    0    |    0    |      0       |
|   B   |    0    |    0    |      L'      |
|   C   |   noe   |    L'   |      0       |
|   D   |   noe   |    L'   |      L'      |
|   E   |   noe   |   2L'   |      L'      |


Bruker eventene A og B:

s_AB^2 = s_AB'^2
dt^2 - dx^2 = dt'^2 - dx'^2
2391[km]^2 - 2678.101[km]^2 = - L'^2
L' = sqrt(2678.101^2[km^2] - 2391^2[km^2]) = 1206.376 [km]




Møtes:


Siden skipet er i ro, midt mellom skip 1 og 2, i romskipenes refferansesystem, \\
er det like langt fra begge skipene til skip 3 og dermed vil hun oppleve det som at lysene ble sendt ut samtidig.

Siden skip 3 beveger seg med samme hastighet som de andre skipene i planetens refferansesystem, \\
vil det si at lysene sees samtidig der også. \\
Dette betyr derfor at lysene må sendes ut på forskjellig tidspunkt i planetens refferansesystem, \\
siden skipene der beveger seg i samme retning som lyset, og laseren fra skip 1 vil bruke lenger tid på å treffe skip 2.

Det er da laseren fra event A som må ha blitt sendt først, i planetens refferansesystem, for at de skal \\
møtes på midten samtidig, siden denne må bevege seg samme vei som skipene, mens laseren fra event B vil \\
gå motsatt vei av skipene, og derfor trenger en veldig mye kortere tid for å nå midtpunktet mellom skipene. \\


"""

c = 3e8
T_B = 7.97*10**(-3)
H_s = 0.891*c
I_s_p2 = 547.72*10**3
#print (I_s_p2 + H_s*T_B + c*T_B)/(H_s + c)
#print - (547.72*10**3)/(H_s - c)
#print H_s*0.0089354856337*10**-3
#print (I_s_p2 + H_s*0.0167498470948)*10**-3
#print (I_s_p2 + H_s*0.00797)*10**-3
#print  



#MittSolSystem.part2A_5(7, friend_seed = 15622)
#MittSolSystem.part2A_8()
