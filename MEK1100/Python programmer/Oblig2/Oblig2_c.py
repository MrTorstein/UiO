#coding=utf-8
from Oblig2_a import *

def plotsqs():
	x1 = [x[0, 34], x[0, 69]]; x2 = [x[0, 69], x[0, 34]]
	x1_ = [x[0, 69], x[1, 69]]; x2_ = [x[0, 34], x[1, 34]]
	plt.plot(\
	x1, [y[159], y[159]], "r", x1_, [y[159], y[169]], "g", \
	x2, [y[169], y[169]], "b", x2_, [y[169], y[159]], "black", \
	x1, [y[84], y[84]], "r", x1_, [y[84], y[99]], "g", \
	x2, [y[99], y[99]], "b", x2_, [y[99], y[84]], "black", \
	x1, [y[49], y[49]], "r", x1_, [y[49], y[59]], "g", \
	x2, [y[59], y[59]], "b", x2_, [y[59], y[49]], "black")

if __name__ == "__main__":
	plt.quiver(x[0: -1: 10, 0: -1: 10],y[0: -1: 10, 0: -1: 10],u[0: -1: 10, 0: -1: 10],v[0: -1: 10, 0: -1: 10])
	plt.plot(xit, yit, "go")
	plotsqs()
	plt.title("Vektorplott for hastighetsfeltet")
	plt.savefig("oblig2_c.png")
	
	plt.show()

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig2_c.py

"""