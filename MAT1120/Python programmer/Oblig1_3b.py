from sympy import Matrix
a = Matrix([[-6, 3, 2, 0], [5, -5, 2, 0], [1, 2, -4, 0]])
print a.rref()

"""
C:\Users\Torstein\Documents\UiO\Mat1120\Python programmer>python Oblig1_3b.py
(Matrix([
[1, 0, -16/15, 0],
[0, 1, -22/15, 0],
[0, 0,      0, 0]]), (0, 1))

"""