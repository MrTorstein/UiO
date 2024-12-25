#coding=utf-8
from Oblig2_4d import *

def oppg4e(a, b, c, d, K, N, alpha, beta):
	from matplotlib.pylab import pcolormesh, zeros, array, zeros_like, show, savefig
	def C(i, j):
		return oppg4d((x[0][i], x[1][j]), (alpha, beta), N)
	x = zeros([2, K])
	l = zeros([K, K])
	for i in xrange(0, K):
		for j in xrange(0, K):
			x[0][i] = (a + i*(b - a))/float(K)
			x[1][i] = (c + j*(d - c))/float(K)
			l[i][j] = C(i, j)
	pcolormesh(l)
	savefig("Oblig2_4e.png")
	show()

if __name__ == "__main__":
	oppg4e(-1, 1, -1, 1, 800, 800, -0.8025, 0.156)