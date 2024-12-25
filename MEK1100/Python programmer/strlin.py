#coding=utf-8
from streamfun import streamfun
from matplotlib.pylab import contour, show, figure, title, savefig

X, Y, psi = streamfun(5)
X_, Y_, psi_ = streamfun(30)

figure()
contour(X, Y, psi)
title("Stomlinjene for (cos(x)sin(y)i - sin(x)cos(y)j), med 5 punkter")
savefig("strlin1.png")
figure()
contour(X_, Y_, psi_)
title("Stomlinjene for (cos(x)sin(y)i - sin(x)cos(y)j, med 30 punkter")
savefig("strlin2.png")
show()