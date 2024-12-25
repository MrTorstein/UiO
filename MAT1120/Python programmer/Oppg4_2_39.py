#coding=utf-8
from numpy import array
from sympy import Matrix
print Matrix(array([[5, 1, 2, 2, 0], [3, 3, 2, -1, -12], [8, 4, 4, -5, 12], [2, 1, 1, 0, -2]])).rref()

#radreduserer matrisen		[   5,   1,   2,   2,   0]
#							[   3,   3,   2,  -1, -12]
#							[   8,   4,   4,  -5,  12]
#							[   2,   1,   1,   0,  -2]