from numpy import linalg, array

I = array([[9, - 3, 0], [- 3, 9, 0], [0, 0, 6]])
print linalg.eig(I)