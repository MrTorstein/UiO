from math import sqrt

xverdier = [1, 1 - sqrt(2)]
for i in range(99):
	"xn+2 = 2x_(n+1) + x_n = 2*1+1 = 3"
	x_n = xverdier[i]
	x_n_1 = xverdier[i+1]
	x_n_2 = 2*x_n_1 + x_n
	print "%3d" % (i+2), x_n_2
	xverdier.append(x_n_2)

"""
Har bare tatt med noen av verdiene fra outputten
Terminalen>oblig1matinf_oppg1b.py
  2 0.171572875254
  3 -0.0710678118655
  4 0.0294372515229
  5 -0.0121933088198
 21 -1.28890269568e-08
 22 -5.21949727883e-09
 23 -2.33280215145e-08
 42 -0.407698098894
 43 -0.984270279703
 44 -2.3762386583
 45 -5.7367475963
 97 -4.60258035988e+20
 98 -1.11116119267e+21
 99 -2.68258042134e+21
100 -6.47632203535e+21
"""