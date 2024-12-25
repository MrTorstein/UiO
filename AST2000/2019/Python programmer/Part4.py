#coding=utf-8
from AST2000SolarSystem import AST2000SolarSystem
from PIL import Image
from numpy import load, pi, array, sqrt, sin, cos, arctan, arcsin, zeros, uint8, sum as Sum, matrix, linspace, meshgrid, save, linalg
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Part1 import part1
seed = 12827
MittSolSystem = AST2000SolarSystem(seed)

def funk(bilde):
	with open(bilde, "rb") as infile:
		himmelkule = load(infile)
	
	alfa_phi = 70*pi/180.; alfa_theta = 70*pi/180.
	theta0 = pi/2.

	x_max = 2*sin(alfa_phi/2.)/(1 + cos(alfa_phi/2.))
	x_min = - x_max
	y_max = 2*sin(alfa_theta/2.)/(1 + cos(alfa_theta/2.))
	y_min = - y_max
	
	dx = 2*x_max/640. + 0.000000001
	dy = 2*y_max/480. + 0.000000001
	
	x = linspace(x_min, x_max, 640)
	y = linspace(y_min, y_max, 480)
	
	X, Y = meshgrid(x, y)
	XY = zeros((480, 640, 2))
	XY[:, :, 0] = X
	XY[:, :, 1] = Y
	
	theta0 = pi/2
	
	liste = zeros((360, 480, 640, 3), dtype = uint8)
	for i in xrange(360 - 1):
		print i
		phi0 = i*pi/180.
		rho = sqrt(X**2 + Y**2)
		c = 2*arctan(rho/2.)
		teller = X*sin(c)
		nevner = rho*sin(theta0)*cos(c) - Y*cos(theta0)*sin(c)
		theta = pi/2. - arcsin(cos(c)*cos(theta0) + teller/rho)
		phi = phi0 + arctan(X*sin(c)/(rho*sin(theta0)*cos(c) - Y*cos(theta0)*sin(c)))
		for n, (j, v) in enumerate(zip(theta, phi)):
			for m, (k, w) in enumerate(zip(j, v)):
				pixnr = AST2000SolarSystem.ang2pix(k, w)
				temp = himmelkule[pixnr]
				liste[i, n, m] = (temp[2], temp[3], temp[4])
	return liste

#save("projeksjon.npy", funk("himmelkule.npy"))


def lsf(filename):
	bilde = Image.open(filename)
	bilde = array(bilde)
	piks = load("projeksjon.npy")
	ki = zeros(360)
	for i in xrange(360):
		ki[i] = Sum((bilde - piks[i])**2)
	print ki[270]
	phi0_best = ki.argmin()
	print "Orienteringsvinkel:", phi0_best, "grader"

#ki[274] = 22330181.0 men denne er riktig =(
#ki[270] = 21156235.0

#lsf("sample0000.png")

def comp_v(dL_sat_1, dL_sat_2):
	get_ref = array(MittSolSystem.get_ref_stars())
	phi_1 = get_ref[0, 0]					#grader
	dL_ref_1 = get_ref[0, 1]				#nm
	phi_2 = get_ref[1, 0]					#grader
	dL_ref_2 = get_ref[1, 1]				#nm
	L0 = 656.3								#nm
	c = 3.e8								#m/s
	
	c = c*(60*60*24*365.24)/149597870700	#Gjør om fra [m/s] til [AU/år]
	phi_1 = phi_1*pi/180					#Gjør om fra grader til radianer
	phi_2 = phi_2*pi/180					#Gjør om fra grader til radianer
	
	dL1 = dL_ref_1 - dL_sat_1
	dL2 = dL_ref_2 - dL_sat_2
	
	""" L - L0)/L0 = vr/c """
	vr1 = (dL1/L0)*c				
	vr2 = (dL2/L0)*c
	
	
	A = matrix([[sin(phi_2), -sin(phi_1)], [-cos(phi_2), cos(phi_1)]])
	vr = matrix([[vr1], [vr2]])
	
	v = (A*vr)/sin(phi_2 - phi_1)
	print "v:", v[0, 0], v[1, 0]


def PosSat(fil = "pos.npy", t = 0.1):
	
	#x**2 + y**2 = s**2
	#y = sqrt(s**2 - x**2)
	"""
	cos(3*pi/4.), sin(3*pi/4)
	r = sqrt(cos(3*pi/4.)**2 + sin(3*pi/4)**2)
	"""

	innfil = open(fil, "rb")
	avstander = load(innfil)
	
	infile = open("planet_positions.npy", "rb")
	ppos, tider = load(infile)
	pf = interp1d(tider, ppos)
	
	p0 = avstander[0]
	p1 = avstander[1]
	p2 = avstander[2]
	p3 = avstander[3]
	p4 = avstander[4]
	p5 = avstander[5]
	p6 = avstander[6]
	p7 = avstander[7]
	s = avstander[8]
	
	
	xp0 =  pf(t)[0, 0]		#AU
	yp0 =  pf(t)[1, 0]		#AU
	xp1 =  pf(t)[0, 1]		#AU
	yp1 =  pf(t)[1, 1]		#AU
	
	
	N = 1.5*10**4			#antall iterasjoner
	tol = 10**-4			#nøyaktighet for posisjonen
	dx0 = -2*s/N			#endring i x retning for en ring rundt stjernen med radius lik avstanden til satelitten
	dx1 = -2*p0/N			#endring i x retning for en ring rundt planet 0 med radius lik avstanden til satelitten
	dx2 = -2*p1/N			#endring i x retning for en ring rundt planet 2 med radius lik avstanden til satelitten

	x0 = s
	x1 = p0
	x2 = p1
	
	num = 0
	x0liste = x0 + dx0*linspace(0, N, N)
	x1liste = x1 + dx1*linspace(0, N, N)
	x2liste = x2 + dx2*linspace(0, N, N)
	y0pliste = sqrt(s**2 - x0liste**2)
	y0mliste = - y0pliste
	y1pliste = sqrt(p0**2 - x1liste**2)
	y1mliste = - y1pliste
	y2pliste = sqrt(p1**2 - x2liste**2)
	y2mliste = - y2pliste

	x1liste = x1liste + xp0
	x2liste = x2liste + xp1
	y1pliste = y1pliste + yp0
	y1mliste = y1mliste + yp0
	y2pliste = y2pliste + yp1
	y2mliste = y2mliste + yp1
	
	x0_1 = 0; y0_1 = 0
	for i in xrange(int(N)):
		x0 = x0liste[i]
		y0p = y0pliste[i]
		y0m = y0mliste[i]
		for j in xrange(int(N)):
			x1 = x1liste[j]
			y1p = y1pliste[j]
			y1m = y1mliste[j]
			if abs(x0 - x1) <= tol:
				if abs(y0p - y1p) <= tol:
					if num == 1:
						x0_2 = x0
						y0_2 = y0m
						num = 2
						break
					x0_1 = x0
					y0_1 = y0p
					num = 1
				elif abs(y0m - y1m) <= tol:
					if num == 1:
						x0_2 = x0
						y0_2 = y0m
						break
					x0_1 = x0
					y0_1 = y0m
					num = 1
		if num == 2:
			break
	if num < 2:
		print "Posisjon til satelitt:", x0_1, y0_1
	tol1 = 9*10**-1
	tol2 = 5*10**-1
	passer = 0
	for i in xrange(int(N)):
		x2 = x2liste[i]
		y2p = y2pliste[i]
		y2m = y2mliste[i]
		
		if passer == 1:
			pass
		
		elif (abs(x2 - x0_1) <= tol1 and abs(y2p - y0_1) <= tol2) or (abs(x2 - x0_1) <= tol1 and abs(y2m - y0_1) <= tol2):
			print "Posisjon til satelitt:", x0_1, y0_1
			passer = 1
		
		elif (abs(x2 - x0_2) <= tol1 and abs(y2p - y0_2) <= tol2) or (abs(x2 - x0_2) <= tol1 and abs(y2m - y0_2) <= tol2):
			print "Posisjon til satelitt:", x0_2, y0_2
			passer = 1
	if passer == 0:
		print "Posisjoner til satelitt (vet ikke hvem):", x0_1, y0_1, x0_2, y0_2


if __name__ == "__main__":
	lsf("find_orient.png")
	comp_v(-0.051748610737, -0.019431558295)
	PosSat("pos.npy")