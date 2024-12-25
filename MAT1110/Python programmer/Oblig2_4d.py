#coding=utf-8
from math import sqrt
from numpy import linalg
from sys import exit

def oppg4d(x_, c, N):
	a = c[0]
	b = c[1]
	M = (1 + sqrt(1 + 4*sqrt(a**2 + b**2)))/2.
	print M
	x = []; x.append([x_[0], x_[1]])
	
	while linalg.norm(x[-1]) <= M:
		x.append([x[-1][0]**2 - x[-1][1]**2 + a, 2*x[-1][0]*x[-1][1] + b])
	
	l = [N+1]
	for j in xrange(len(x)):
			if linalg.norm(x[j]) >= M:
				l.append(j)
	n = min((N, min(l)))
	return n/float(N)

if __name__ == "__main__":
	n = oppg4d((0.325, 0.35), (-0.8, 0.156), 800)
	print "%.4f" %(n)

"""
C:\Users\Torstein\Documents\UiO\Mat1110\Python programmer>Oblig2_4d.py
0.0138

"""