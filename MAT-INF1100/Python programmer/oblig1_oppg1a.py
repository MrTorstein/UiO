xverdier = [1, 1]
for i in range(99):
	"xn+2 = 2x_(n+1) + x_n = 2*1+1 = 3"
	x_n = xverdier[i]
	x_n_1 = xverdier[i+1]
	x_n_2 = 2*x_n_1 + x_n
	print "%3d" % (i+2), x_n_2
	xverdier.append(x_n_2)

"""
Har bare tatt med noen av resultatene
Terminalen>oblig1matinf_oppg1a.py
  2 3
  3 7
  4 17
  5 41
 97 6733044458057842709277507685523012161
 98 16255007246704249599863612909970798723
 99 39243058951466341909004733505464609607
100 94741125149636933417873079920900017937
"""