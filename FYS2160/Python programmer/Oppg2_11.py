"""Funker ikke"""
from numpy import linspace, exp
from scipy import misc
from matplotlib.pyplot import figure, plot, legend, title, xlabel, ylabel, show

N = 100						# Antall partikler
q = 80						# Total energimengde
d = 1000					# Antall steg

A_q = linspace(0, q, d)
A_omg = (A_q + N - 1) ** (A_q + N - 1) / ( A_q ** A_q * ( N - 1 ) ** ( N - 1 ) )


figure()

plot( A_q, A_omg / 10**49 )
plot( A_q, linspace(1, 2, d) )

show()

