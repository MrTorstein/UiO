with open("running.txt", "r") as infile:
	t = []
	v = []
	for line in infile:
		tnext, vnext = line.strip().split(',')
		t.append(float(tnext))
		v.append(float(vnext))

from matplotlib.pyplot import subplot, plot, show, xlim, ylim, \
title, savefig
from numpy import zeros

a = zeros(len(v))
s = zeros(len(v))
for i in xrange(len(v) - 1):
	a[i] = (v[i + 1] - v[i])/(t[i + 1] - t[i])
	s[i + 1] = s[i] + (v[i + 1] + v[i])/2.*(t[i + 1] - t[i])

subplot(2, 1, 1)
plot(t , a)
xlim(-100, 7000)
title("Akselerasjon")

subplot(2, 1, 2)
plot(t, s)
title("Strekning")

savefig("oblig2oppg1c.png")
show()

"""
C:\Users\Torstein\Documents\UiO\Inf1100\Python programmer>oblig2matinf_oppg1c.py

"""