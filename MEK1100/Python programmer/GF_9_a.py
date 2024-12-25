#coding=utf-8
from numpy import linspace, zeros, meshgrid
from matplotlib.pyplot import show, figure, quiver
from mpl_toolkits.mplot3d import Axes3D

x = linspace(-20, 20, 4)
y = linspace(-20, 20, 4)
z = 1200 + 250./4*x*y
"""
zx = 1200 + 250./4*x*1
zy = 1200 + 250./4*1*y
a = [zx, zy]
"""
"""
z = meshgrid(x, y)
quiver(x, y, z)
show()
"""
"""
fig = figure()
ax = fig.add_subplot(111, projection = "3d")
ax.contour(x, y, a)
show()
"""

fig = figure()
ax = fig.add_subplot(111, projection = "3d")
ax.quiver(x, y, z, x, y, z)
show()