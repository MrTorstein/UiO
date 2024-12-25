#coding=utf-8
from matplotlib.pylab import polyfit, polyval, linspace, array, plot, show, scatter, xlabel, ylabel, savefig

class MakeFunk():
	
	def __init__(self, x, y):
		self.x = x
		self.y = y
	
	def Eval(self, ny_x, deg):
		self.ny_x = ny_x
		p = polyfit(self.x, self.y, deg)
		y_res = polyval(p, ny_x)
		
		return p, y_res
	
	def Plot(self):
		plot(self.ny_x, y_res)
		scatter(x, y)
		xlabel("x")
		ylabel("y")
		savefig("Halleffekt.png")
		show()

if __name__ == "__main__":
	x = array([18, 66, 114, 164, 218, 272])
	y = array([1.8, 6.4, 10.8, 22.6, 31.8, 38.7])
	Funkmaker = MakeFunk(x, y)
	poly, y_res = Funkmaker.Eval(x, 1)
	Funkmaker.Plot()
	#print poly
	#0.15373869
	
	q = 1.6*10**-19
	RH = 0.15373869*1*10**-3/(25*10**-3)
	print "RH =", RH
	N = 1/(RH*q)
	print "N =", N