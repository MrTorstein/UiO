from sympy import *
a = Matrix([[1, 1, 1, 100], \
[1, (-1/2.-sqrt(3)*sqrt(1)/2.), (-1/2.+sqrt(3)*sqrt(1)/2.), 0], \
[1,(-1/2.+sqrt(3)*sqrt(1)/2.) , (-1/2.-sqrt(3)*sqrt(1)/2.), 0]]).rref()
for i in range(len(a)- 1):
	for j in range(len(a[i])):
		a[i][j] = float(a[i][j])

print a

"""
C:\Users\Torstein\Documents\UiO\Mat1110\Python programmer>Radreduksjon.py
(Matrix([
[1.0, 0.0, 0.0, 33.3333333333333],
[0.0, 1.0, 0.0, 33.3333333333333],
[0.0, 0.0, 1.0, 33.3333333333333]]), [0, 1, 2])

"""