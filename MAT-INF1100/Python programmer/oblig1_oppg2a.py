def binom(n, i):
	b = 1
	for j in range(1, n-i+1):
		n = float(n)
		i = float(i)
		a = (i+j)/j
		b = b*a
	return b

print binom(9998, 4)
print binom(100000, 70)
print binom(1000, 500)

"""
Terminalen>oblig1matinf_oppg2a.py
4.16083629103e+14
8.14900007814e+249
2.70288240945e+299
"""