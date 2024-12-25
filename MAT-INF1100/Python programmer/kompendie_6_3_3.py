def diff_2(g, f_0, f_1, f_2, a0, a1, a2, N):
	x_ppp = a0
	x_pp = a1
	x_p = a2
	for i in range(2, N+1):
		x = g(i - 2) + f_2(i - 2)*x_ppp + f_0(i - 2)*x_pp + f_1(i - 2)*x_p
		print x
		x_ppp = x_pp
		x_pp = x_p
		x_p = x

def g(n):
	return 0
def f_0(n):
	return 1
def f_1(n):
	return 1
def f_2(n):
	return 1

diff_2(g, f_0, f_1, f_2, 0, 1, 1, 10)

"""
C:\Users\Torstein\Documents\UiO\Inf1100\Python programmer>
kompendie_oppg6_3_3.py
2
4
7
13
24
44
81
149
274
"""