#coding=utf-8
from numpy import *
from matplotlib.pyplot import * 

#
def force(dx,d,k,q):
    if dx<d:
        F = k*abs(dx-d)**q
    else:
        F = 0.0
    return F

N = 3    # nr of balls
m = 1. # kg
k = 20000 # N/m
q = 1.
d =  1. # m
v0 = 1. # m/s
time = 0.025 # s
dt = 0.00001 # s

n = int(round(time/dt))
x = zeros((n,N),float)
v = x.copy()
t = zeros(n,float)
# Initial conditions

for i in range(N):
    x[0,i] = d*i

v[0,0] = v0

for i in range(n-1):
    # Find force in vector F
    F = zeros(N,float)
    for j in range(1,N):
        dx = x[i,j] - x[i,j-1]
        F[j] = F[j] + force(dx,d,k,q)
    for j in range(N-1):
        dx = x[i,j+1] - x[i,j]
        F[j] = F[j] - force(dx,d,k,q)
    # Euler-Cromer vectorized step
    a = F/m
    v[i+1] = v[i] + a*dt
    x[i+1] = x[i] + v[i+1]*dt
    t[i+1] = t[i] + dt

for j in range(N):
    plot(t,v[:,j])
    if j==0:
        hold('on')
    if j==N-1:
        hold('off')

print 'v/v0 = ',v[n-1,:]/v0
#savefig("oblig5_k.png")
show()