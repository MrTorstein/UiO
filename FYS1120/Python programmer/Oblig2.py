#coding=utf-8
import time
start = time.clock()
import sys
from matplotlib.pylab import pi, cos, array, linspace, zeros, cross, plot, show, title, legend, xlabel, ylabel, savefig, \
figure, sqrt, axis
class Sim():
	"""
	felt = en matrise med E og/eller B i streng først og feltverdiene i sine respektive plasser
	st = starttid
	et = sluttid(endtime)
	sp = startposisjon
	sh = starthastighet
	dt = integralsteg
	fortegn = fortegn på ladning
	"""
	def __init__(self, felt):
		if len(felt[0]) == 2:
			if felt[0][0] == "E":
				if callable(felt[1][0]):
					self.E = felt[1][0]
				else:
					self.E = lambda t: felt[1][0]
				self.B = felt[1][1]
			elif felt[0][0] == "B":
				if callable(felt[1][0]):
					self.E = felt[1][1]
				else:
					self.E = lambda t: felt[1][1]
				self.B = felt[1][0]
			else:
				print "Noe er feil med hvordan feltet er oppgitt \n", felt, "\n Skulle vært på formen [('E', 'B'), ((1, 1, 1), (0, 0, 1))]"
				sys.exit()
		
		elif len(felt[0]) == 1:
			if felt[0] == "E":
				if callable(felt[1]):
					self.E = felt[1]
				else:
					self.E = lambda t: felt[1]
				self.B = 0
			elif felt[0] == "B":
				self.E = lambda t: 0
				self.B = felt[1]
			else:
				print "Noe er feil med hvordan feltet er oppgitt \n", felt, "\n Skulle vært på formen ['E' eller 'B', (1, 1, 1)]"
				sys.exit()
	
	def Initialverdier(self, st = 0,  sp = (0, 0, 0), sh = (0, 0, 0), fortegn = -1, m = 9.11*10**(-31), d = 0, r = 0):
		try:
			abs(fortegn) == 1
		except:
			fortegn = int(fortegn/float(fortegn))
		self.st = st										#s
		self.sp = sp										#m
		self.sh = sh										#m/s
		self.q = fortegn*1.60*10**-19						#C
		self.m = m											#kg
		self.d = d
		self.r = r
	
	def a(self, t, p, h, d = 0, r = 0):
		E = self.E
		B = self.B
		q = self.q
		if d != 0:
			if abs(p) <= 0.5*d:
				FE = self.q*E(t)
			else:
				FE = 0
		else:
			FE = self.q*E(t)
		
		if r != 0:
			if abs(p) <= r:
				FB = self.q*cross(h, B)
			else:
				FB = 0
		else:
			FB = q*cross(h, B)
		
		return (FE + FB)/self.m
	
	def Simuler(self, et, dt):
		#tol = 10**-10
		st = self.st
		a = self.a
		N = int((et - st)/dt)
		t = linspace(st, et, N)
		h = zeros((N, 3))
		h[0] = self.sh
		p = zeros((N, 3))
		p[0] = self.sp
		
		bredde = 20											#bredden paa lastingsmaaleren
		sys.stdout.write("Laster... [%s]" % (" "*bredde))	#
		sys.stdout.flush()
		sys.stdout.write("\b" * (bredde + 1)) # returnerer til starten av linja, etter "["
		
		mellomrom = len(t)/bredde
		teller = 0
		
		for i in xrange(N - 1):
			h[i + 1] = h[i] + self.a(t[i], p[i, 0], h[i], self.d, self.r)*dt
			p[i + 1] = p[i] + h[i + 1]*dt
			#if abs(p[i + 1, 0]) <= tol and abs(p[i + 1, 1]) <= tol and abs(p[i + 1, 2]) <= tol:
			#	print t[i]									#1.78871924795e-11 sek
			if i/mellomrom > teller:
				sys.stdout.write("%s" %("="))
				teller = i/mellomrom
			
		print "]"
		print "Done!"
		return h, p, t


from mpl_toolkits.mplot3d import Axes3D
from numpy.ma import masked_inside, masked_outside


B = ("B", array((0, 0, 2)))

S = Sim(B)
S.Initialverdier(sh = (10**4, 0, 0))
h1, p1, t1 = S.Simuler(30*10**(-12), 10**(-15))

figure()
plot(t1, p1)
legend(["x(t)", "y(t)", "z(t)"], loc = "best")
xlabel("t[s]")
ylabel("x[m]")
title("Posisjoner")
#savefig("2a_Posisjoner_2D.png")

figure()
plot(t1, h1)
legend(["$v_x(t)$", "$v_y(t)$", "$v_z(t)$"], loc = "best")
xlabel("t[s]")
ylabel("v[$m/s$]")
title("Hastigheter")
#savefig("2a_Hastigheter_2D.png")

fig = figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(p1[:, 0], p1[:, 1], p1[:, 2])
ax.set_xlabel("Posisjon [m]")
ax.set_ylabel("Posisjon [m]")
ax.set_zlabel("Posisjon [m]")
title("Banen til partikkelen")
#savefig("2a_Posisjoner_3D.png")

show()

#print pi*9.11*10**(-31)/(1.60*10**-19)							#1.78874431714e-11 sek

S.Initialverdier(sh = (5*10**3, 0, 0.2*10**3))
h2, p2, t2 = S.Simuler(30*10**(-12), 10**(-15))

fig = figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(p2[:, 0], p2[:, 1], p2[:, 2])
ax.set_xlabel("Posisjon [m]")
ax.set_ylabel("Posisjon [m]")
ax.set_zlabel("Posisjon [m]")
title("Banen til partikkelen")
#savefig("2d_Posisjoner_3D.png")

show()


E = lambda t: array((25*10**3/(90*10**-6)*cos(1.602177*10**-19*1.69/(1.67262158*10**-27)*t), 0, 0))
Felt = [("E", "B"), (E, array((0, 0, 1.69)))]

S = Sim(Felt)

S.Initialverdier(fortegn = 1, m = 1.67262158*10**-27, d = 9.10**-5, r = 0)
h, p, t = S.Simuler(et = 3*10**-7, dt = 10**-13)

figure()
plot(p[:, 0], p[:, 1])
legend(["Bane til partikkel"], loc = "best")
xlabel("x[m]")
ylabel("y[m]")
title("Syklotronsimulasjon")
#savefig("3a_SyklotronSim.png")
show()


S.Initialverdier(fortegn = 1, m = 1.67262158*10**-27, d = 9.10**-5, r = 50*10**-3)
h, p, t = S.Simuler(et = 1.5*10**(-6), dt = 10**(-12))

figure()
plot(t, p)
legend(["x(t)", "y(t)", "z(t)"], loc = "best")
xlabel("t[s]")
ylabel("x[m]")
title("Posisjoner Syklotronsimulasjon")
savefig("3b_SyklotronSim_pos.png")

figure()
plot(t, h)
legend(["$v_x(t)$", "$v_y(t)$", "$v_z(t)$"], loc = "best")
xlabel("t[s]")
ylabel("x[m]")
title("Hastigheter Syklotronsimulasjon")
savefig("3b_SyklotronSim_hast.png")
show()

#print h[-1]													#[526268.01870168, -8094112.2943176, 0]
#print sqrt(h[-1, 0]**2 + h[-1, 1]**2 + h[-1, 2]**2)			#8111202.86151 m/s


S.Initialverdier(fortegn = 1, m = 1.67262158*10**-27, d = 9.10**-5, r = 1)
h, p, t = S.Simuler(et = 10**-5, dt = 10**-10)

figure()
plot(t, p)
legend(["x(t)", "y(t)", "z(t)"], loc = "best")
xlabel("t[s]")
ylabel("x[m]")
title("Posisjoner Syklotronsimulasjon")
axis([0, 3.5*10**-6, -1.5, 1.5])
#savefig("3e_SyklotronSim_pos.png")
show()

#print 0.5*1.67262158*10**-27*(sqrt(h[-1, 0]**2 + h[-1, 1]**2 + h[-1, 2]**2))**2 #2.19285616009e-11 J

end = time.clock()
print "Tiden programmet brukte er:", end - start, "sek", "eller", (end - start)/60., "min"
