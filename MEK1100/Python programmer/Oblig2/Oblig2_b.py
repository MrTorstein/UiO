#coding=utf-8
from Oblig2_a import *

plt.figure(1)
plt.contourf(x, y, u, v, levels = np.linspace(-500, 500, 200))
plt.colorbar()
plt.plot(xit, yit, "go")
plt.title("Hastighetsfelt fra -500 til 500")
plt.savefig("oblig2_b1.png")

plt.figure(2)
plt.plot(xit, yit, "go")
plt.contourf(x, y, u, v, levels = np.linspace(500, 4500, 200))
plt.colorbar()
plt.title("Hastighetsfelt fra 500 til 4500")
plt.savefig("oblig2_b2.png")

plt.show()

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig2_b.py

"""