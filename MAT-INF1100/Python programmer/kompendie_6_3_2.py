def diff_2(g, f_0, f_1, a0, a1, N):
	x_pp = a0
	x_p = a1
	for i in range(2, N+1):
		x = g(i - 2) + f_0(i - 2)*x_pp + f_1(i - 2)*x_p
		print x
		x_pp = x_p
		x_p = x
def g(n):
	return 0
def f_0(n):
	return 1
def f_1(n):
	return 1

diff_2(g, f_0, f_1, 0, 1, 10)

"""
C:\Users\Torstein\Documents\UiO\Inf1100\Python programmer>
kompendie_oppg6_3_2.py
1
2
3
5
8
13
21
34
55
"""