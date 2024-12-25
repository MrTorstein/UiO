#coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

nIterasjoner = 100
N = 40
V = np.zeros((N, N))

"""Grenseverdier for potensialet"""
linfun = np.linspace(0, 1, N)
topp = 1 - np.cos(linfun * (3 * np.pi + np.pi / 2))
V[0, :] = 0													# Unødvendig linje
V[N - 1, :] = 10
V[:, 0] = 10 * linfun
V[:, N - 1] = 10 * topp

fig = plt.figure(1)
ax = fig.add_subplot(111, projection = '3d') 				# Tilgang til 3d-plots
plt.ion() 													# Få se stegene i den iterative løseren
vx, vy = np.meshgrid(range(0, N), range(0, N))

"""Iterasjoner for å regne ut potensialet i området"""
for it in xrange(0, nIterasjoner):
	for ix in xrange(1, N - 1):
		for iy in xrange(1, N - 1):
			V[iy, ix] = 1. / 4 * (V[iy + 1, ix] + V[iy, ix + 1] + V[iy - 1, ix] + V[iy, ix - 1])
	plt.cla() 												# Fjerne gammelt plot så ikke figuren begynner å henge
	ax.plot_surface(vx, vy, V, cmap = cm.viridis, shade = False)
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('V(x,y)')
	ax.set_title("Iterasjon %d/%d" % (it + 1, nIterasjoner))
	plt.pause(0.001)


Ey, Ex = np.gradient(- V)
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.quiver(vx, vy, Ex, Ey, units = 'width', scale = 100) 	# scale er tilpasset for at plottet skal se fint ut
ax2.set_title("Elektrisk felt")
plt.show(block=True)
