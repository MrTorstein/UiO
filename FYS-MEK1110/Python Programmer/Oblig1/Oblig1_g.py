def Euler(A):
	from numpy import zeros
	t = zeros(10000)
	a = zeros(10000); a[0]= A(0, 0)
	v = zeros(10000)
	s = zeros(10000)
	s_ = s[0]
	t_ = t[0]
	v_ = v[0] 
	i = 0
	a_ = 1
	dt = 1
	
	while a_ >= 1E-13:
		a_ = A(t_, v[i])
		t_ += dt
		s_ = s[i] + v_*(dt)
		v_ = v[i] + a_*(dt)
		a[i + 1] = a_
		v[i + 1] = v_
		s[i + 1] = s_
		t[i + 1] = t_
		i += 1
	return t, s, v, a

if __name__ == "__main__":	#ikke saa mye betydning
	from Oblig1_e import A	#importerer funksjonen A fra forrige oppgave
	t, s, v, a = Euler(A)
	v = v[v > 0]
	print "Toppfarten til loperen er %.0fm/s" %v[-1]

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>Oblig1_g.py
Toppfarten til loperen er 34m/s

"""
