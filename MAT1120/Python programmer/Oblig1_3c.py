#coding=utf-8
def Rekkeutvikler(v, matrise, antall):	#Definerer en funksjon for aa regne ut utviklinga
	"""
	v		= Matrise med utgangspunktet for fordelinga
	matrise	= Utviklingsmatrisa P
	antall	= Liste over de "dagene" vi er ute etter 	
	"""
	slutt = max(antall)			#Finner den siste dagen vi er ute etter
	res = []
	for i in xrange(slutt):		#Kjorer utregninga
		v = matrise*v
		if i in antall or i == slutt - 1:
			res.append(v)
	return res

from numpy import matrix, array
P = matrix([[0.4, 0.3, 0.2], [0.5, 0.5, 0.2], [0.1, 0.2, 0.6]])
xyz_0 = matrix([[20000], [25000], [8000]])
dager = [4, 10, 20]

resultater = Rekkeutvikler(xyz_0, P, dager)
for i in xrange(len(dager)):	#Printer
	print "Dag %2d:" %dager[i],"\n", resultater[i], "\n"

"""
C:\Users\Torstein\Documents\UiO\Mat1120\Python programmer>python Oblig1_3c.py
Dag  4:
[[ 16033.59]
 [ 22077.34]
 [ 14889.07]]

Dag 10:
[[ 16000.21314027]
 [ 22000.49081422]
 [ 14999.29604551]]

Dag 20:
[[ 16000.00010775]
 [ 22000.00024812]
 [ 14999.99964414]]

"""