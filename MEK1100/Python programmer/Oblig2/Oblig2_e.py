#coding=utf-8
from Oblig2_c import *

dudy = np.gradient(u, axis = 0)
dvdx = np.gradient(v, axis = 1)
curlz = dvdx - dudy 

if __name__ == "__main__":
	plt.contourf(x, y, curlz, levels = np.linspace(-1000, 1000, 13))
	plt.colorbar()
	plt.plot(xit, yit, "go")
	plotsqs
	plt.title("Virvlinga til u$i$ + v$j$")
	plt.savefig("oblig2_e.png")
	
	plt.show()

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig2_e.py

"""