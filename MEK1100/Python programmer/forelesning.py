#coding=utf-8
from numpy import *
from matplotlib.pyplot import *

t = linspace(-5, 5, 11)
x, y = meshgrid(t, t)
quiver(x, y, x)
show()