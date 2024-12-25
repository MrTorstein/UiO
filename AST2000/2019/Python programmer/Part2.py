#coding=utf-8
import sys
from numpy import sqrt, zeros, zeros_like, random, sum as Sum, pi, exp, linspace, round, array, save, load, linalg
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, savefig
from AST2000SolarSystem import AST2000SolarSystem
seed = 12827
MittSolSystem = AST2000SolarSystem(seed)

N = int(2e6)
tid = 82.

MTS = MittSolSystem.star_mass
MTP = MittSolSystem.mass
AP = MittSolSystem.number_of_planets

PTP = zeros((AP, 2, N))
HTP = zeros_like(PTP)
temp1 = array(zip(MittSolSystem.x0, MittSolSystem.y0))
temp2 = array(zip(MittSolSystem.vx0, MittSolSystem.vy0))

PTP[:, :, 0] = temp1[:, :]
HTP[:, :, 0] = temp2[:, :]

dt = tid/N
G = 4*pi**2

print "Regner ut planetbaner"
bredde = 20												#bredden på lastingsmåleren
sys.stdout.write("Laster... [%s]" % (" "*(bredde - 1)))
sys.stdout.flush()
sys.stdout.write("\b" * (bredde)) 						#returnerer til starten av linja, etter "["
mellomrom = (N)/bredde
teller = 0


a_x_0 = - G*MTS*PTP[:, 0, 0]/sqrt(PTP[:, 0, 0]**2 + PTP[:, 1, 0]**2)**3
a_y_0 = - G*MTS*PTP[:, 1, 0]/sqrt(PTP[:, 0, 0]**2 + PTP[:, 1, 0]**2)**3
HTP[:, :, 0] = HTP[:, :, 0] + array(zip(a_x_0, a_y_0))*dt/2. - 0.00001


for i in xrange(N - 1):
	a_x = - (G*MTS*PTP[:, 0, i])/sqrt(PTP[:, 0, i]**2 + PTP[:, 1, i]**2)**3
	a_y = - (G*MTS*PTP[:, 1, i])/sqrt(PTP[:, 0, i]**2 + PTP[:, 1, i]**2)**3
	HTP[:, :, i + 1] = HTP[:, :, i] + array(zip(a_x, a_y))*dt
	PTP[:, :, i + 1] = PTP[:, :, i] + HTP[:, :, i + 1]*dt
	
	if i/mellomrom > teller:
		sys.stdout.write("%s" %("="))
		teller = i/mellomrom
			
print "\nDone!"



plot(\
PTP[0, 0, ::1000], PTP[0, 1, ::1000], \
PTP[1, 0, ::1000], PTP[1, 1, ::1000], \
PTP[2, 0, ::1000], PTP[2, 1, ::1000], \
PTP[3, 0, ::1000], PTP[3, 1, ::1000], \
PTP[4, 0, ::1000], PTP[4, 1, ::1000], \
PTP[5, 0, ::1000], PTP[5, 1, ::1000], \
PTP[6, 0, ::1000], PTP[6, 1, ::1000], \
PTP[7, 0, ::1000], PTP[7, 1, ::1000])

legend(["Planet 1", "Planet 2", "Planet 3", "Planet 4", "Planet 5", "Planet 6", "Planet 7", "Planet 8"], loc = "best")
savefig("Planetbaner.png")
xlabel("[AU]")
ylabel("[AU]")
show()

"""
t = linspace(0, tid, N)
E = zeros((AP, N - 1))
M = MTS + MTP
mu = MTS*MTP/(M)
for i in xrange(N - 1):
	rA = sqrt((PTP[:, 0, i])**2 + (PTP[:, 1, i])**2)
	rB = sqrt((PTP[:, 0, i + 1])**2 + (PTP[:, 1, i + 1])**2)
	dr = rB - rA
	dt = t[i + 1] - t[i]
	v = dr/dt
	E[:, i] = 1/2.*mu*v**2 - G*M*mu/rB
print "Total energi i starten:\n", E[:, 0]
print "Total energi til slutt:\n", E[:, -1]
print "Diff: \n", E[:, 0] - E[:, -1]
"""

xy = zeros((2, AP, N))
for i in xrange(AP):
	for j in xrange(2):
		xy[j, i, :] = PTP[i, j, :]
PTP = xy

MittSolSystem.orbit_xml(PTP[:, :, ::1000], linspace(0, tid, N)[::1000])

MittSolSystem.check_planet_positions(PTP, tid, N/tid)
#Får plutselig ikke laget .npy fila lenger, pga. MemoryError, men har allerede gjordt det