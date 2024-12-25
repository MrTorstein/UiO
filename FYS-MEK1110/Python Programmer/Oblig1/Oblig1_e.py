""" Definerer en funksjon for aa finne akselrasjonen """
def A(t, v):
	F = 400	#N
	rho = 1.293	#kg/m^3
	A = 0.45	#m^2
	CD = 1.2
	m = 80.	#kg
	return (F - 0.5*rho*A*CD*v**2)/m	#m/s^2

""" Definerer en funksjon for aa bruker Eulers metode """
def Euler(A):
	t = [0]
	a = [A(0, 0)]
	v = [0]
	s = [0]
	s_ = s[0]
	t_ = t[0]
	i = 0
	dt = 0.01
	
	while s_ < 100:
		a_ = A(t_, v[i])
		t_ += dt
		v_ = v[i] + a_*(dt)
		s_ = s[i] + v_*(dt)
		a.append(a_)
		v.append(v_)
		s.append(s_)
		t.append(t_)
		i += 1
	return t, s, v, a

if __name__ == "__main__":	#bare for at jeg skal kunne bruke funksjoner videre
	from matplotlib.pyplot import plot, legend, show, xlabel, savefig
	t, s, v, a = Euler(A)
	
	plot(t, s, t, v, t, a)
	legend(["s(t)", "v(t)", "a(t)"], loc = "best")
	xlabel("t")
	savefig("Oblig1_e.png")
	show()


"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>Oblig1_e.py

"""