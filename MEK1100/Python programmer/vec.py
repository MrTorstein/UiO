#coding=utf-8
from velfield import velfield
from matplotlib.pylab import quiver, show, title, savefig

x, y, u, v = velfield(21)
quiver(x, y, u, v)
title("Vektorplott for (cos(x)sin(y)i - sin(x)cos(y)j), med 21 punkter")
savefig("vec.png")
show()