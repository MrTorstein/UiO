#coding=utf-8
import time
import sys
start = time.clock()
from matplotlib.pylab import cos, linspace, zeros, cross, plot, show, title, legend, xlabel, ylabel, savefig, figure, sqrt, axis
def Sim(sp = (0, 0, 0), sh = (0, 0, 0), st = 0, et = 3*10**-7, dt = 10**-13, d = 9.10**-5, r = 0):
	"""
	st = starttid
	et = sluttid(endtime)
	sp = startposisjon
	sh = starthastighet
	dt = integralsteg
	d = vally gap i syklotron
	r = unnslipningsradius til syklotron
	"""
	
	q = 1.602177*10**-19									#C
	m = 1.67262158*10**-27 									#kg	
	N = int((et - st)/dt)									#lengde paa integrasjon
	t = linspace(st, et, N)									#sek
	h = zeros((N, 3))										#m/s
	h[0] = sh
	p = zeros((N, 3))										#m
	p[0] = sp
	w = 1.602177*10**-19*1.69/(1.67262158*10**-27)			#sek^-1
	B = (0, 0, 1.69)										#T
	
	bredde = 15
	sys.stdout.write("Laster... [%s]" % (" "*bredde))
	sys.stdout.flush()
	sys.stdout.write("\b" * (bredde + 1)) # return to start of line, after '['
	
	mellomrom = len(t)/bredde
	teller = 0
	for i in xrange(N - 1):
		if abs(p[i, 0]) <= 0.5*d:
			FE = lambda t: (25*10**3/(90*10**-6)*cos(w*t)*q, 0, 0)
		else:
			FE = lambda t: 0
		FB = q*cross(h[i], B)
		if abs(p[i, 0]) <= r and r != 0:
			FB = q*cross(h[i], B)
		elif r != 0:
			FB = 0
		else:
			FB = q*cross(h[i], B)
		a = (FE(t[i]) + FB)/m
		h[i + 1] = h[i] + a*dt
		p[i + 1] = p[i] + h[i + 1]*dt
		if i/mellomrom > teller:
			sys.stdout.write("%s" %("="))
			teller = i/mellomrom
			
	print "]"
	print "Done!"
	return h, p, t
"""
h, p, t = Sim()

figure()
plot(p[:, 0], p[:, 1])
legend(["Bane til partikkel"], loc = "best")
xlabel("x[m]")
ylabel("y[m]")
title("Syklotronsimulasjon")
#savefig("3a_SyklotronSim.png")

h, p, t = Sim(r = 50*10**-3)
figure()
plot(t, p)
legend(["x(t)", "y(t)", "z(t)"], loc = "best")
xlabel("t[s]")
ylabel("x[m]")
title("Posisjoner Syklotronsimulasjon")
#savefig("3b_SyklotronSim_pos.png")

figure()
plot(t, h)
legend(["$v_x(t)$", "$v_y(t)$", "$v_z(t)$"], loc = "best")
xlabel("t[s]")
ylabel("x[m]")
title("Hastigheter Syklotronsimulasjon")
#savefig("3b_SyklotronSim_hast.png")

h, p, t = Sim(et = 1.45*10**(-6), dt = 10**(-12), r = 50*10**-3)
figure()
plot(p[:, 0], p[:, 1])
legend(["Bane til partikkel"], loc = "best")
xlabel("x[m]")
ylabel("y[m]")
title("Syklotronsimulasjon")
show()

#print h[-1]													#[526268.01870168, -8094112.2943176, 0]
#print sqrt(h[-1, 0]**2 + h[-1, 1]**2 + h[-1, 2]**2)			#8111202.86151 m/s
"""
h, p, t = Sim(et = 7.8*10**-6, dt = 3*10**-11, r = 1)

figure()
plot(t, p)
legend(["x(t)", "y(t)", "z(t)"], loc = "best")
xlabel("t[s]")
ylabel("x[m]")
title("Posisjoner Syklotronsimulasjon")
axis([0, 8*10**-6, -1.5, 1.5])
#savefig("3e_SyklotronSim_pos.png")
show()

print 0.5*1.67262158*10**-27*(sqrt(h[-1, 0]**2 + h[-1, 1]**2 + h[-1, 2]**2))**2

end = time.clock()
print "Tiden programmet brukte er:", end - start, "sek", "eller", (end - start)/60., "min"