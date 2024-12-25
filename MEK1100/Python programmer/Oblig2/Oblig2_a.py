#coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import h5py


# Dette ser ut til å virke med versjon 7.3 MAT-filer
with h5py.File("data73.mat", "r") as file:
	X = file["x"]
	Y = file["y"]
	U = file["u"]
	V = file["v"]
	Xit = file["xit"]
	Yit = file["yit"]
	x = np.transpose(X[:,:])
	y = np.transpose(Y[:,:])
	u = np.transpose(U[:,:])
	v = np.transpose(V[:,:])
	xit = np.transpose(Xit[:])
	yit = np.transpose(Yit[:])

if __name__ == "__main__":
	print x.shape, y.shape, u.shape, v.shape, xit.shape, yit.shape, "lengde i (y, x) retning"
	print x[0][1] - x[0][0], y[1][0] - y[0][0], y[0][0], y[-1][0], "avstanden mellom punktene og start og sluttpunkt for y aksen"

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig2_a.py
(201L, 194L) (201L, 194L) (201L, 194L) (201L, 194L) (1L, 194L) (1L, 194L) lengde i (y, x) retning
0.5 0.5 -50.0 50.0 avstanden mellom punktene og start og sluttpunkt for y aksen

"""