from numpy import linspace, zeros_like, pi
from matplotlib.pyplot import plot, show, xlabel, ylabel, title, legend, figure, savefig

"""
Likning for dempet fjerpendel:
z''(t) = - b / m * z'(t) - k / m * z(t)


Da har vi disse variabler aa leke med:
b		= Dempningsfaktor
k		= Fjerkonstant
m		= Masse til objekt festet paa fjera
z_0		= Startposisjon for pendel
v_0		= Starthastighet for pendel
T		= Sluttid
N		= Antall tidssteg
"""

def RK4(k, l, m, z_0, v_0, T, N):
	
	"""
	a(t) = k * v(t) + l * z(t) + m * t
	
	Variablene blir da:
	 - k
	 - l
	 - m
	 - N
	 - T
	
	Likningssettene for Runge Kutta 4 blir:
	x_1,n = z[n]
	v_1,n = v[n]
	a_1,n = - b / m * v_1,n - k / m * x_1,n
	
	x_2,n = x_1,n + v_1,n * dt / 2
	v_2,n = v_1,n + a_1,n * dt / 2
	a_2,n = - b / m * v_2,n - k / m * x_2,n
	
	x_3,n = x_1,n + v_2,n * dt / 2
	v_3,n = v_1,n + a_2,n * dt / 2
	a_3,n = - b / m * v_3,n - k / m * x_3,n
	
	x_4,n = x_1,n + v_3,n * dt
	v_4,n = v_1,n + a_3,n * dt
	a_4,n = - b / m * v_4,n - k / m * x_4,n
	
	vsnitt_n = 1 / 6 * ( v_1,n + v_2,n + v_3,n + v_4,n )
	asnitt_n = 1 / 6 * ( a_1,n + a_2,n + a_3,n + a_4,n )
	
	z[n + 1] = z[n] + vsnitt_n * dt
	v[n + 1] = v[n] + asnitt_n * dt
	"""
	
	dt = T / N
	def f(v_, x_, t_):
		return k * v_ + l * x_ + m * t_
	
	t = linspace(0, T, N)
	z = zeros_like(t); z[0] = z_0
	v = zeros_like(t); v[0] = v_0
	
	for i in xrange(len(t) - 1):
		vsn = 0
		asn = 0
		t_ = t[i] + dt / 2.
		
		x_1 = z[i]
		v_1 = v[i]
		a_1 = f(v_1, x_1, t[i])
		
		x_2 = x_1 + v_1 * dt / 2.
		v_2 = v_1 + a_1 * dt / 2.
		a_2 = f(v_2, x_2, t_)
		
		x_3 = x_1 + v_2 * dt / 2.
		v_3 = v_1 + a_2 * dt / 2.
		a_3 = f(v_3, x_3, t_)
		
		x_4 = x_1 + v_3 * dt
		v_4 = v_1 + a_3 * dt
		a_4 = f(v_4, x_4, t[i + 1])
		
		vsn = (v_1 + 2*v_2 + 2*v_3 + v_4) / 6.
		asn = (a_1 + 2*a_2 + 2*a_3 + a_4) / 6.
		
		z[i + 1] = z[i] + vsn * dt
		v[i + 1] = v[i] + asn * dt
	
	return t, z, v

def plotte(var, funk, N, legends = 0):
	for i in xrange(len(var)):
		figure(i)
		plot(var[i], funk[i])
		title("Funksjon med %d antall tidssteg" %(N * 10**(i)))
		xlabel("tid [s]")
		ylabel("posisjon [m]")
		if legends != 0:
			legend(["%s" %(legends[i])])
		#savefig("Figur_%04d" %(i))
	show()
	
if __name__ == "__main__":
	
	b 	= 0.1	#kg/s
	k 	= 10	#N/m
	m 	= 0.1	#kg
	z_0	= 0.1	#m
	v_0	= 0		#m/s
	T	= 10	#s
	N	= 1e2 
	
	t1, z1, v1 = RK4(- b / m, - k / m, 0, z_0, v_0, T, N)
	t2, z2, v2 = RK4(- b / m, - k / m, 0, z_0, v_0, T, N * 1e1)
	t3, z3, v3 = RK4(- b / m, - k / m, 0, z_0, v_0, T, N * 1e2)
	t4, z4, v4 = RK4(- b / m, - k / m, 0, z_0, v_0, T, N * 1e3)


	plotte([t1, t2, t3, t4], [z1, z2, z3, z4], N)
	
	dt = T / (N * 1e2)
	i = pi/dt
	print int(i), z3[int(i)]
	
"""
Kjore eksempel:
C:\Users\Torstein\Documents\UiO\Fys2130\Python programmer>python Oblig2_4_4.py
3141 0.0207258161169

"""