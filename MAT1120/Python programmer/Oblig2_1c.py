from sympy import Matrix
a = Matrix([[3, 1, 1, 5, 3], [1, 2, 0, 3, 2], [1, 0, 1, 2, 1], [5, 3, 2, 10, 6]])
print a.rref()