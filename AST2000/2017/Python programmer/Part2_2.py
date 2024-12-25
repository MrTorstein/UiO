#coding=utf-8
from AST2000SolarSystem import AST2000SolarSystem
from numpy import sqrt, zeros, zeros_like, array, pi, linspace, random, cos, sum as Sum, where
from matplotlib.pyplot import plot, legend, show, figure, savefig, xlabel, ylabel, title
import sys
seed = 12827
MittSolSystem = AST2000SolarSystem(seed)
random.seed(seed)

N = 10**5
tid = 105
"""
liste = [("Nr.", "    Masse     ", "      Avstand      ")]
liste2 = [sqrt(MittSolSystem.x0[0]**2 + MittSolSystem.y0[0]**2), sqrt(MittSolSystem.x0[1]**2 + MittSolSystem.y0[1]**2), sqrt(MittSolSystem.x0[2]**2 + MittSolSystem.y0[2]**2), sqrt(MittSolSystem.x0[3]**2 + MittSolSystem.y0[3]**2), sqrt(MittSolSystem.x0[4]**2 + MittSolSystem.y0[4]**2), sqrt(MittSolSystem.x0[5]**2 + MittSolSystem.y0[5]**2), sqrt(MittSolSystem.x0[6]**2 + MittSolSystem.y0[6]**2), sqrt(MittSolSystem.x0[7]**2 + MittSolSystem.y0[7]**2)]
for i in xrange(8):
	laveste = 0
	for j in xrange(8):
		if liste2[j] < liste2[laveste]:
			if j != laveste:
				laveste = j
	liste.append((laveste, MittSolSystem.mass[laveste], liste2[laveste]))
	liste2[laveste] = 100
for i in liste:
	print i
"""
#Bruker planet nr. 2

MTS = MittSolSystem.star_mass
MTP = MittSolSystem.mass[2]

PTO = zeros((3, 2, N))
HTO = zeros((2, 2, N))
PCM = sqrt(MittSolSystem.x0[2]**2 + MittSolSystem.y0[2]**2)*MTP/(MTP + MTS)

temp1 = array((sqrt(MittSolSystem.x0[2]**2 + MittSolSystem.y0[2]**2), 0))
temp2 = array((0, sqrt(MittSolSystem.vx0[2]**2 + MittSolSystem.vy0[2]**2)))
temp3 = array((PCM, 0))

PTO[0, :, 0] = temp3
PTO[2, :, 0] = temp1
HTO[1, :, 0] = temp2

dt = tid/float(N)
G = 4*pi**2

a_x_S_0 = - G*MTP*(PTO[1, 0, 0] - PTO[2, 0, 0])/sqrt((PTO[1, 0, 0] - PTO[2, 0, 0])**2 + (PTO[1, 1, 0] - PTO[2, 1, 0])**2)**3
a_x_P_0 = - G*MTS*(PTO[2, 0, 0] - PTO[1, 0, 0])/sqrt((PTO[2, 0, 0] - PTO[1, 0, 0])**2 + (PTO[2, 1, 0] - PTO[1, 1, 0])**2)**3
a_y_S_0 = - G*MTP*(PTO[1, 1, 0] - PTO[2, 1, 0])/sqrt((PTO[1, 0, 0] - PTO[2, 0, 0])**2 + (PTO[1, 1, 0] - PTO[2, 1, 0])**2)**3
a_y_P_0 = - G*MTS*(PTO[2, 1, 0] - PTO[1, 1, 0])/sqrt((PTO[2, 0, 0] - PTO[1, 0, 0])**2 + (PTO[2, 1, 0] - PTO[1, 1, 0])**2)**3
HTO[:, :, 0] = HTO[:, :, 0] + array(zip((a_x_S_0, a_x_P_0), (a_y_S_0, a_y_P_0)))*dt/2.

bredde = 20												#bredden på lastingsmåleren
sys.stdout.write("Laster... [%s]" % (" "*(bredde - 1)))
sys.stdout.flush()
sys.stdout.write("\b" * (bredde)) 						#returnerer til starten av linja, etter "["
mellomrom = (N)/bredde
teller = 0

for i in xrange(0, N - 1):
	a_x_S = - G*MTP*(PTO[1, 0, i] - PTO[2, 0, i])/sqrt((PTO[1, 0, i] - PTO[2, 0, i])**2 + (PTO[1, 1, i] - PTO[2, 1, i])**2)**3
	a_x_P = - G*MTS*(PTO[2, 0, i] - PTO[1, 0, i])/sqrt((PTO[2, 0, i] - PTO[1, 0, i])**2 + (PTO[2, 1, i] - PTO[1, 1, i])**2)**3
	a_y_S = - G*MTP*(PTO[1, 1, i] - PTO[2, 1, i])/sqrt((PTO[1, 0, i] - PTO[2, 0, i])**2 + (PTO[1, 1, i] - PTO[2, 1, i])**2)**3
	a_y_P = - G*MTS*(PTO[2, 1, i] - PTO[1, 1, i])/sqrt((PTO[2, 0, i] - PTO[1, 0, i])**2 + (PTO[2, 1, i] - PTO[1, 1, i])**2)**3
	HTO[:, :, i + 1] = HTO[:, :, i] + array(zip((a_x_S, a_x_P), (a_y_S, a_y_P)))*dt
	PTO[1 : 3, :, i + 1] = PTO[1 : 3, :, i] + HTO[:, :, i + 1]*dt
	PTO[0, :, i + 1] = array(((MTS*PTO[1, 0, i + 1] + MTP*PTO[2, 0, i + 1])/(MTP + MTS), MTS*(PTO[1, 1, i + 1] + MTP*PTO[2, 1, i + 1])/(MTP + MTS)))
	
	if i/mellomrom > teller:
		sys.stdout.write("=")
		teller = i/mellomrom

print ("]  Done!")
figure()
plot(PTO[0, 0, ::100], PTO[0, 1, ::100], PTO[1, 0, ::100], PTO[1, 1, ::100], PTO[2, 0, ::100], PTO[2, 1, ::100])
legend(["CM", "Stjerne", "Planet"], loc = "best")
xlabel("[AU]")
ylabel("[AU]")
savefig("Baner.png")

t = linspace(0, tid, N)
figure()
plot(t[::100], HTO[0, 1, ::100]*1e3)
title("Hastighet til Stella i x-retning")
xlabel("tid [aar]")
ylabel("hastighet [AU/aar] 10^-3")
savefig("Solhastighet.png")

Sy = HTO[0, 1]
Sy = Sy + random.normal(0, 1/5.*Sy.max(), N)
figure()
plot(t[::100], Sy[::100]*1e3)
title("Hastighet til Stella i x-retning, med gaussisk forsyrrelse")
xlabel("tid [aar]")
ylabel("hastighet [AU/aar] 10^-3")
savefig("Solhastighet_forstyrring.png")

show()

E = zeros(N)
mu = MTS*MTP/(MTS + MTP)
for i in xrange(N - 1):
	rA = sqrt(((PTO[2, 0, i] - PTO[1, 0, i])**2 + (PTO[2, 1, i] - PTO[1, 1, i])**2))
	rB = sqrt(((PTO[2, 0, i + 1] - PTO[1, 0, i + 1])**2 + (PTO[2, 1, i + 1] - PTO[1, 1, i + 1])**2))
	dr = rB - rA
	dt = t[i + 1] - t[i]
	v = dr/dt
	E[i] = 1/2.*mu*v**2 - G*(MTP + MTS)*mu/rB

print "Total energi til system ved ti forskjellige tidspunkter, gjevnt fordelt\n", E[::N/10]


S = array(zip(t, Sy))
#plot(S[::100, 0], S[::100, 1])
#show()


def dsum(t0, P, vr):
	vmod = vr*cos(2*pi/P*(S[:, 0] - t0))
	summen = Sum((S[:, 1] - vmod)**2)
	return t0, P, vr, summen
liste = []

for t0_ in xrange(20, 30):
	for P_ in xrange(45, 66):
		for vr_ in xrange(50, 71):
			liste.append(dsum(t0_, P_, vr_*1e-6))

liste = array(liste)
#print liste[where(liste[:, -1] == min(liste[:, -1]))]
# 22.000000 60.000000 0.000050 0.001065
t0 = 22.
P = 60.
vr = 0.00005
t = linspace(0, tid, N)
vmod = vr*cos(2*pi/P*(t - t0))

mp = MTS**(2/3.)*vr*P**(1/3.)/(2*pi*G)**(1/3.)
print mp # 5.70561302559e-05
print "Differansen =", mp - MTP
"""
C:\Users\Torstein\Documents\UiO\Ast2000\Python programmer>python Part2_2.py
Differansen = 2.25665134698e-05

"""