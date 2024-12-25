#coding=utf-8
from Oblig2_c import *

dudx = np.gradient(u, axis = 1)
dvdy = np.gradient(v, axis = 0)
div = dudx + dvdy

if __name__ == "__main__":
	plt.contourf(x, y, div, levels = np.linspace(-1000, 1000, 13))
	plt.colorbar()
	plt.plot(xit, yit, "go")
	plotsqs()
	plt.title("Divergensen til u$i$ + v$j$")
	plt.savefig("oblig2_d.png")
	
	plt.show()

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig2_d.py

"""