#coding=utf-8
from sympy import *
a = Matrix([\
[-2, 1, 0, 0, 1, 0], [1, -2, 1, 0, 0, 0], [0, 1, -2, 1, 0, 0], [0, 0, 1, -2, 1, 0], [1, 0, 0, 1, -2, 0]]).rref()
for i in range(len(a)- 1):
	for j in range(len(a[i])):
		a[i][j] = float(a[i][j])

print a

"""
C:\Users\Torstein\Documents\UiO\Mat1110\Python programmer>Oblig2_2a.py
(Matrix([
[1.0, 0.0, 0.0, 0.0, -1.0, 0.0],
[0.0, 1.0, 0.0, 0.0, -1.0, 0.0],
[0.0, 0.0, 1.0, 0.0, -1.0, 0.0],
[0.0, 0.0, 0.0, 1.0, -1.0, 0.0],
[0.0, 0.0, 0.0, 0.0,  0.0, 0.0]]), [0, 1, 2, 3])

"""