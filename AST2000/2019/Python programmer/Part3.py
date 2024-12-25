#coding=utf-8
from AST2000SolarSystem import AST2000SolarSystem
seed = 12827
MittSolSystem = AST2000SolarSystem(seed)
from numpy import array, linspace, zeros, zeros_like, sqrt, arccos, linalg, sum as Sum, pi, exp, load, append, dot
from matplotlib.pyplot import plot, show, axis, legend, xlabel, ylabel, figure, savefig
from scipy.interpolate import interp1d
from sys import exit

class Solcelle():
	
	def __init__(self, Pn = 8):
		self.Pn = Pn
	
		self.ST = MittSolSystem.temperature
		self.SR = MittSolSystem.star_radius*10**3					#Stjerneradius						m
		self.PR = MittSolSystem.radius[Pn - 1]*10**3				#Radius til planet					m
		self.XP = sqrt(MittSolSystem.x0[Pn - 1]**2 + MittSolSystem.y0[Pn - 1]**2)
		self.XP = self.XP*149597870700 - self.SR					#Avstand fra stjerne til planet		m
		
		self.N = 100000
		self.lam = linspace(10**-9, 3000*10**-9, self.N)			#Bølgelengder
		self.h = 6.63*10**(-34)										#Plancks konstant					J s
		self.k = 1.38064852*10**(-23)								#Boltzmanns konstant				m^2 kg s^-2 K^-1
		self.c = 3*10**8											#Lyshastigheten						m/s
		self.sig = 5.67*10**-8										#Stefan-Boltzmanns konstant			W m^-2 K^-4
		self.PT_Ol_min = 260										#Mintemp for at liv på planet skal overleve				K
		self.PT_Ol_max = 390										#Maxtemp for at liv på planet skal overleve				K
	
	def Wien(self):
		self.lam_max = 2.9*10**(-3)/self.ST						#Wiens forskyvningslov: lambda_max = k/T
		print self.lam_max, "m"									#3.15352555705e-07	m
	
	def B(self):
		B = (2*self.h*self.c**2)/self.lam**5*1/(exp((self.h*self.c)/(self.k*self.ST*self.lam)) - 1)
		plot(self.lam*10**9, B, [315, 315.3], [0, 3*10**14])
		xlabel("m")
		show()
	
	def RegnUt(self, printes = []):
		""" """													#	Verdier for planet nr.8
		FS = self.sig*self.ST**4								#	405498853.415			W/m^2
		L = FS*4*pi*self.SR**2									#	9.0758828188*10^27		W
		FP = L/(4*pi*self.XP**2)								#	521.971150301			W/m^2
		A = 40/(FP*0.12)										#	0.638604898262			m^2
		E = FP*2*pi*self.PR**2									#	1.25641626297e*10^16	W
		PT = (E/(4*pi*self.PR**2*self.sig))**(1/4.)				#	260.470390623			K
		
	    
		X_Ol_max = sqrt(L/(8*pi*self.PT_Ol_min**4*self.sig))	#	1.18055516266*10^12		m
		X_Ol_min = sqrt(L/(8*pi*self.PT_Ol_max**4*self.sig))	#	524691183403.0			m
		
		if "avstander" in printes:
			print "Avstander fra planetene til stjernen:"
			for i in xrange(8):
				avstand = sqrt(MittSolSystem.x0[i]**2 + MittSolSystem.y0[i]**2)
				print "%2i %10e m" %(i + 1, avstand*149597870700 - self.SR)
				#Bare en av planetene faller innenfor "Habitable zone", planet nr. 8
			printes.pop(printes.index("avstander"))
			print " "
		
		if len(printes) > 0:
			"""
			Printbare ting er:
			FS			#Fluks til stjerne									W/m^2
			L           #Luminositet til stjerne							W
			FP          #Fluksen som treffer/sendes ut fra planeten			W/m^2
			A           #Arealet av solcellepanelet							m^2
			E           #Energien, per sekund, som treffer planeten			W
			PT          #Temperaturen til planeten på grunn av stjerna		K
			X_Ol_max	#Maximumavstand for en planet med liv				m
			X_Ol_min    #Minimumavstand for en planet med liv				m
			avstander	#Avstand til planetene								m
			"""
			for i in printes:
				print i, eval(i)





class Simuler:
	
	def __init__(self, starttid = 0, startpos = 0, starthast = [0, 0], drivstoff = 10000, dest = 8):
		self.satelittmasse = 1100								#										[kg]
		self.drivstoff = drivstoff								#										[kg]
		self.starttid = starttid								#Starttid for simulasjon hvis den er gitt, ellers 0
		self.reisetid = 0										#Reisetiden fra utskytningspunkt til destinasjon
		self.sluttid = 82.										#Sluttid for simulasjon, når hjemplanet har beveget seg 20 runder
		self.startpos2 = startpos								#Startposisjon hvis den er gitt			[m]
		self.starthast = starthast								#Starthastighet hvis gitt				[m/s]
		self.dp = 4.4e-9*7e12									#Endring i bevegelsesmengde per tidsenhet
		self.G = 4*pi**2										#Gravitasjonskonstanten 				[AU^3 yr^-2 SM^-1]
		self.N = 1*10**4										#Antall totale steg
		self.dest = dest - 1									#Nummer på destinasjonsplanet
		innfil = open("planet_positions.npy", "rb")
		self.planetpos, tider = load(innfil)
		self.posfunk = interp1d(tider, self.planetpos)			#Funksjon for å finne posisjon til planeter
		self.M = self.satelittmasse + self.drivstoff			#Total masse til raketten				[kg]
		self.km_to_AU = 6.68458712*10**-9						#										[AU/m]
		self.tol = 0.01											#tolleransen jeg skal treffe innenfor	[AU]
	
	def plot_planet(self):
		pp = self.planetpos
		fig1 = figure(1)
		plot(\
		pp[0, 0, ::1000], pp[1, 0, ::1000],\
		pp[0, 1, ::1000], pp[1, 1, ::1000],\
		pp[0, 2, ::1000], pp[1, 2, ::1000],\
		pp[0, 3, ::1000], pp[1, 3, ::1000],\
		pp[0, 4, ::1000], pp[1, 4, ::1000],\
		pp[0, 5, ::1000], pp[1, 5, ::1000],\
		pp[0, 6, ::1000], pp[1, 6, ::1000],\
		pp[0, 7, ::1000], pp[1, 7, ::1000])
		legend(["Planet 1", "Planet 2", "Planet 3", "Planet 4", "Planet 5", "Planet 6", "Planet 7", "Planet 8"], loc = "best")
	
	def Show(self):
		show()
	
	def finnstart(self):
		print "finner start"
		r = MittSolSystem.radius[0]
		pf = self.posfunk
		dest = self.dest
		self.starttid = 0
		
		for i in xrange(int(self.sluttid + 1)):
			p0i = pf(i)[:, 0]
			pdi = pf(i)[:, dest]
			if linalg.norm(p0i - pdi) <= 4 and self.starttid == 0:
				self.starttid = i
				self.startpos = p0i + r*self.km_to_AU*(pdi - p0i)/linalg.norm(pdi - p0i)
				self.reisetid = linalg.norm(pdi - p0i)/1.49481497191			#hentet v fra launchhastigheten
				pdiT = pf(self.starttid + self.reisetid)[:, dest]
				self.startpos2 = p0i + r*self.km_to_AU*(pdiT - p0i)/linalg.norm(pdiT - p0i)
				plot(p0i[0], p0i[1], "mo")										#hjemplanet ved utskytning
				plot(pdi[0], pdi[1], "go")										#destplant ved utskytning
				#plot(pdiT[0], pdiT[1], "yo")									#destplanet ved tid T=reisetid
				"""																#plotter posisjonene til satelitten på overflaten til planeten
				plot(self.startpos[0], self.startpos[1], "rx")
				plot(self.startpos2[0], self.startpos2[1], "bx")
				"""
		else:
			pass
	
	def launch(self):
		self.finnstart()
		pot_energy = self.G*MittSolSystem.mass[0]/linalg.norm(array([3.50934777912 - MittSolSystem.x0[0], 0.00016954538814]))
		v_esc = sqrt(pot_energy*2)
		
		r = MittSolSystem.radius[0]
		pf = self.posfunk
		dest = self.dest - 1
		pdiT = pf(self.starttid + self.reisetid)[:, dest]
		p0i = pf(self.starttid)[:, 0]
		self.retning = (pdiT - p0i)/linalg.norm(pdiT - p0i)
		
		self.starthast = v_esc*self.retning
		self.startpos2 += linalg.norm(array([3.50934777912 - MittSolSystem.x0[0], 0.00016954538814]))*self.retning
		plot(self.startpos2[0], self.startpos2[1], "bx")
		self.starttid += 1200./(60*60*24*365)
		print "Skutt opp"
	
	def boost(self, retning, hast):
		dv = 10*retning
		self.M -= 0.9106026455567189
		return dv

	def simuler(self):
		if type(self.startpos2) == type(0):
			self.launch()
		
		N = self.N; t1 = self.starttid; t2 = self.sluttid; pf = self.posfunk; dest = self.dest; G = self.G
		dt = (t2 - t1)/N
		MO = append(MittSolSystem.star_mass, MittSolSystem.mass)
		PTS = zeros((2, N))
		PTS[:, 0] = self.startpos2
		HTS = self.starthast
		HTS_ = array([0, 0])
		t = linspace(t1, t2, N)
		
		for i in xrange(N - 1):
			p = pf(t[i])[:, dest]
			p_x = append(PTS[0, i], PTS[0, i] - pf(t[i])[0])
			p_y = append(PTS[1, i], PTS[1, i] - pf(t[i])[1])
			a_x = - Sum(G*MO*p_x/sqrt(p_x**2 + p_y**2)**3)
			a_y = - Sum(G*MO*p_y/sqrt(p_x**2 + p_y**2)**3)
			HTS_ = HTS + array([a_x, a_y])*dt
			
			if self.reisetid != 0:
				rt = self.reisetid
				
				if linalg.norm(PTS[:, i]) >= linalg.norm(pf(0)[:, 0]):
					ubrukelig = 1
					plotbar = pf(t[i])[:, dest]
				
				elif linalg.norm(PTS[:, i] + HTS_*(rt - t[i]) - pf(rt)[:, dest]) >= self.tol:
					retning = (PTS[:, i] + HTS_*(rt - t[i]) - pf(t[i] + rt)[:, dest])/linalg.norm(PTS[:, i] + HTS_*(rt - t[i]) - pf(t[i] + rt)[:, dest])
					"""
					if retning[1] < 0:
						retning = array([- retning[1], retning[0]])
					
					else:
						retning = array([retning[1], - retning[0]])
					"""
					plot((PTS[:, i] + HTS_*(rt - t[i]))[0],(PTS[:, i] + HTS_*(rt - t[i]))[1] , "mx")
					plot(pf(t[i] + rt)[0, dest], pf(t[i] + rt)[1, dest], "gx")
					HTS_ += self.boost(retning, HTS_)
				
				elif linalg.norm(PTS[:, i] + HTS_*(rt - t[i]) - pf(rt)[:, dest]) < self.tol:
					ubrukelig = 1
					plotbar = pf(t[i])[:, dest]
				
			else:
				w = arccos(dot(pf(1)[:, dest], pf(0)[:, dest])/(linalg.norm(pf(1)[:, dest])*linalg.norm(pf(0)[:, dest])))	#[1/year]
				P = PTS[:, i]
				nummer = 0
				while linalg.norm(P) < linalg.norm(pf(0)[:, dest]):
					P += HTS_*dt
					nummer += 1
				self.reisetid = nummer*dt
				vink = arccos(dot(p, PTS[:, i] + P)/(linalg.norm(p)*linalg.norm(PTS[:, i] + P)))
				self.reisetid = w/vink
				retning = PTS[:, i]/linalg.norm(PTS[:, i])
				HTS_ += self.boost(retning, HTS_)
			
			PTS[:, i + 1] = PTS[:, i] + HTS_*dt
		plot(plotbar[0], plotbar[1], "ro")
		plot(PTS[0, ::10], PTS[1, ::10], "rx")
		savefig("satelittbane.png")
		
tilfelle = Simuler()
tilfelle.plot_planet()
tilfelle.simuler()
tilfelle.Show()
#av en eller anne grunn er det nå noe galt med planet_positions.npy fila mi, som gjør at programmet ikke fungerer lenger =(