#coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import sys
import random
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from numpy import sqrt, zeros, zeros_like, random, sum as Sum, pi, exp, linspace, sin, cos
from AST2000SolarSystem import AST2000SolarSystem
from ast2000solarsystem_27_v2 import AST2000SolarSystem as AST2000SolarSystem1
from numpy.ma import masked_inside, masked_outside
seed = 12827
MittSolSystem = AST2000SolarSystem(seed)
random.seed(seed)

class part1:
	def __init__(self):
		self.H2gPrMol = 2.01588								# Hydrogen gass weight [g/mol]
		self.N_a = 6.02214179e23							# Avogadros number

		self.mass_sat = 1100								# Satellite/Rocket mass [kg]
		self.k = 1.380648e-23								# Boltzmann's constant [J/K]
		self.mass1Particle = self.H2gPrMol/self.N_a*1e-3	# [kg]

	def oppdatere(self, posisjoner, hastigheter, L):		# Funksjon for å oppdatere partikler som kolliderer med vegger
		m = self.mass1Particle
		a = masked_outside(posisjoner, 0, L)
		hastigheter += -2*hastigheter*a.mask
		
		b = masked_inside(posisjoner[:, 0], L/10., 9*L/10.)
		c = masked_inside(posisjoner[:, 1], L/10., 9*L/10.)
		d = masked_inside(posisjoner[:, 2], 0, L)
		e = masked_outside(posisjoner[:, 2], 0, L)
		hullteller = Sum(e.mask*b.mask*c.mask)
		p = abs(Sum(m*hastigheter[:, 0]*b.mask*c.mask*e.mask))
		return hastigheter, p, hullteller
	
	def simulate_box(self, n, N, dt, box_dim, T):
		"""
		Regner ut bevegelsen til pratiklene	som blir brukt til drivstoff.
		
		Utregningsvariabler:
		n		=		Antall partikler
		N		=		Antall tidssteg
		dt		=		Lengden på ett tidssteg
		box_dim	=		Dimensjonen til boksen
		T		=		Temperaturen til gassen
		"""

		m = self.mass1Particle
		sigma = sqrt(T*self.k/m)						# Standaravvik gaussisk hastighetsfordeling
		mean = 0										# Gjennomsnitsverdi for gausisk fordeling

		pos = random.random((N, 3))*box_dim				# Posisjonsvektorer
		vel = random.normal(mean, sigma, (N, 3))		# Hatighetsvektorer

		x_momentum = 0.0								# Bevegelsesmengde i x-retning
		no_escapes = 0									# Antall partikler som unslipper motoren
		
		print "Simulerer boks"
		bredde = 20												# Bredden på lastingsmåleren
		sys.stdout.write("Laster... [%s]" % (" "*(bredde - 1)))
		sys.stdout.flush()
		sys.stdout.write("\b" * (bredde)) 						# Returnerer til starten av linja, etter "["
		mellomrom = (N)/bredde
		teller = 0
		
		for i in xrange(N):
			pos += vel*dt
			vel, p, tel = self.oppdatere(pos, vel, box_dim)
			no_escapes += tel
			x_momentum += p
			
			if i/mellomrom > teller:
				sys.stdout.write("%s" %("="))
				teller = i/mellomrom
		sys.stdout.write("] Done!\n")
		return x_momentum, no_escapes

	def launch_sim(self, no_boxes, v_end, A, no_particles_sec):
		"""
		Simulerer utskytning

		Variabler:
		no_boxes:           antall drivstoffbokser
		v_end:              den ønskede farten i forhold til overflaten ved slutten av simulasjonen
		A:                  dp/dt for boksen beskrevet i første funksjon
		no_particles_sec:   antall partikler som forlater en boks per sekund

		returnerer:    Hvor mye tid det tar å oppnå ønsket fart og hvor mye drivstoff som er brukt
		"""
		v = 0																# Utgangshastighet til rakett i forhold til overfalten
		time = 0               												# Starttid
		T = 20.*60                                      					# Total tid, 20 min
		Nt = 1200															# Antall tidssteg
		dt = float(T)/Nt                                					# Lengde på tidssteg
		G = 6.67428*10**(-11)												# Newtons gravitasjonskonstant
		
		initial_fuel_mass = 50000                      						# Drivstoffmasse 
		Masse = MittSolSystem.mass[0]*1.989*10**30
		fuel_needed = 0
		counter = 0
		while v < v_end:
			B = self.mass1Particle*no_particles_sec*no_boxes				# dm/dt for boksen i forrige funksjon
			M = self.mass_sat + initial_fuel_mass							# Total masse
			b = A*no_boxes
			seconds = 0
			pos = MittSolSystem.radius[0]*1e3
			fuel_left = initial_fuel_mass
			v = 0
			for i in xrange(Nt):
				dv = - G*Masse/pos**2
				a = b/M + dv
				M -= dt*B
				fuel_left -= dt*B
				if fuel_left <= 0:
					print "Tom for drivstoff"
					break
				v = v + a*dt
				if v <= 0:
					break
					
				pos += v*dt
				v_end = sqrt(2*G*Masse/pos)
				seconds += 1

			if v > v_end:
				print "Naadde maalet"
				return v_end, no_boxes, pos, seconds, initial_fuel_mass, initial_fuel_mass - fuel_left
				
			if fuel_left <= 0:
				print "legge til mer drivstoff"
				initial_fuel_mass *= 1.001
			
			counter += 1
			if counter > 5000:
				print "Simulation will never reach Escape Velocity"
				break
			no_boxes = no_boxes*1.001

	def boost_sim(self, no_boxes,v_start, v_end, A, no_particles_sec, fuel_mass):
		"""
		Simulates a boost

		Variables:
		no_boxes:           the number of fuel boxes you need
		v_start:			th initial velocity
		v_end:              the desired velocity at the end of launch
		A:                  dp/dt for the box described in first method
		no_particles_sec:   number of particles pr second leaving the box in first method
		fuel_mass:			fuelmass at start
		
		returns:    How much time it takes to achieve the boost and how much fuel you have used up.
		"""
		
		A = A*no_boxes				                                  		# dp/dt for the whole engine
		B = self.mass1Particle*no_particles_sec								# dm/dt for the box described in first method
		v = v_start															# initial velocity relative to surface of planet	
		time = 0               												# initialize time variable
		sluttid = 1000.														# sek
		Nt = 100000															# Number of time steps
		dt = sluttid/Nt                                						# Time step

		initial_fuel_mass = fuel_mass                      					# Calculate/set fuel mass
		M = self.mass_sat + initial_fuel_mass           					# Total mass
		fuel_needed = 0
		for i in xrange(Nt):
			""" Here you need to update the new value of the velocity aswell as the new total mass value M """ 
			v += (A/M)*dt
			M -= B*no_boxes*dt
			time += dt
			fuel_needed = initial_fuel_mass - (M - self.mass_sat)
			if M < self.mass_sat:
				print "tom for drivstoff"
				print M - self.mass_sat
				print time, B*no_boxes*dt
				print v, v_end
				return 0, 0
			elif v >= v_end:
				""" Boost is successful. Save values and end method """
				#fuel_needed = ...
				print "Boost Succesful"
				return fuel_needed, time

if __name__ == "__main__":
	tilfelle = part1()
	
	N_ = 2.3e4
	t = 1e-10
	p, anpart = tilfelle.simulate_box(n = 10**5, N = int(N_), dt = t/N_, box_dim = 10**(-6), T = 10**4)
	print "p = ", p					#1.19705406766e-19
	print "anpart =", anpart		#7607
	anpart = anpart/t
	p = p/t
	print p
	print anpart
	
	
	R = MittSolSystem.radius[0]*1e3
	M = MittSolSystem.mass[0]*1.989*10**30
	G = 6.67428*10**(-11)
	v_esc = sqrt((2*G*M)/R)

	anpart = 7.607e13
	p = 1.19705406766e-9
	
	
	v_end, n_box, pos, sec, start_fuel, fuel_needed = tilfelle.launch_sim(no_boxes = 3e+13, v_end = v_esc, A = p, no_particles_sec = anpart)
	
	print tilfelle.boost_sim(no_boxes = n_box, v_start = 0, v_end = 1000, A = p, no_particles_sec = anpart, fuel_mass = 1000)
	
	"""
	drivs_trengs = 403.07425431874253
	tid_brukt = 3.0199999999999796
	"""
	
	init_sat_pos = (MittSolSystem.x0[0] + MittSolSystem.radius[0]*6.68458712*10**-9, MittSolSystem.y0[0])
	MittSolSystem.engine_settings(p, n_box, anpart, start_fuel, sec, init_sat_pos, 0)
	leggtily = MittSolSystem.vy0[0]*sec/(60*60*24*365) + MittSolSystem.radius[0]*6.68458712*10**-12*sin(2*pi/MittSolSystem.period[0]*sec/(60*60*24))
	calculated_y = init_sat_pos[1] + leggtily
	leggtilx = v_end*(sec/4.)*6.68458712*10**-12 - (MittSolSystem.radius[0] - MittSolSystem.radius[0]*cos(2*pi/MittSolSystem.period[0]*sec/(60*60*24)))*6.68458712*10**-12
	calculated_x = init_sat_pos[0] + leggtilx
	print "x,y", calculated_x, calculated_y						# 3.50933553617 6.2191647674e-05
	pos_after_launch = calculated_x, calculated_y
	a = MittSolSystem.mass_needed_launch(pos_after_launch, test = True)
