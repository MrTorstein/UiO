#a)
from sympy import pretty
from numpy import matrix, array, linalg, pi, sin, cos
def MKM(datasett, funk):
	temp1 = []
	temp2 = []
	for i in xrange(len(datasett)):
		temp1.append(funk(datasett[i][0]))
		temp2.append([datasett[i][1]])
	
	X = matrix(temp1)
	y = matrix(temp2)
	
	X_T = X.transpose()
	
	X_TX = X_T*X
	X_Ty = X_T*y
	
	B = linalg.inv(X_TX)*X_Ty
	b = ()
	for i in xrange(len(B)):
		b = b + (B[i, 0],)
	return b

datasett = [[0.08, 4.05], [0.12, 4.15], [0.20, 3.85], [0.38, -0.22]]
def funk1(data):
	return [data**0, data**1, data**2]
B1 = MKM(datasett, funk1)
print "funk: y = %.2f + %.2fx + %.2fx^2" % B1



#b)
def funk2(data):
	return [sin(2*pi*data), cos(2*pi*data)]
B2 = MKM(datasett, funk2)
print "funk: y = %.2fsin(2*pi*x) + %.2fcos(2*pi*x)" % B2



#c)
x = array(array(datasett)[:, 0])
y = array(array(datasett)[:, 1])

def kurve1(B, x):
	return array(B[0]*x**0 + B[1]*x**1 + B[2]*x**2)

epsilon1 = y - kurve1(B1, x)
norm1 = linalg.norm(epsilon1)


def kurve2(B, x):
	return array(B[0]*sin(2*pi*x) + B[1]*cos(2*pi*x))

epsilon2 = y - kurve2(B2, x)
norm2 = linalg.norm(epsilon2)

print "||epsilon|| for kurve 1 =", norm1
print "||epsilon|| for kurve 2 =", norm2