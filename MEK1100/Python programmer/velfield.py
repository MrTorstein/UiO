#coding=utf-8
from numpy import linspace, sin, cos, meshgrid
def velfield(n):
	x = linspace(-5, 5, n)
	y = linspace(-5, 5, n)
	X, Y = meshgrid(x, y)
	vx = cos(X)*sin(Y)
	vy = -sin(X)*cos(Y)
	return X, Y, vx, vy

