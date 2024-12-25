def A(t, v):
	F = 400		#N
	fv = 25.8	#sN/m
	m = 80		#kg
	return (F - fv*v)/m	#m/s^2

if __name__ == "__main__":
	from Oblig1_g import Euler	#henter Euler fra forrige oppgave
	t, s, v, a = Euler(A)
	v = v[v > 0]
	print "Toppfarten til loperen er %.0fm/s" %v[-1]

"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>Oblig1_h.py
Toppfarten til loperen er 16m/s

"""