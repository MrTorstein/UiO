#coding=utf-8
from Oblig2_d import *

r1 = - sum(v[159, 34:70])*0.5; g1 = sum(u[159:169, 69])*0.5; b1 = + sum(v[169, 34:70])*0.5; s1 = - sum(u[159:169, 34])*0.5
r2 = - sum(v[ 84, 34:70])*0.5; g2 = sum(u[ 84: 99, 69])*0.5; b2 = + sum(v[ 99, 34:70])*0.5; s2 = - sum(u[ 84: 99, 34])*0.5
r3 = - sum(v[ 49, 34:70])*0.5; g3 = sum(u[ 49: 59, 69])*0.5; b3 = + sum(v[ 59, 34:70])*0.5; s3 = - sum(u[ 49: 59, 34])*0.5

sirkkur1 = r1 + g1 + b1 + s1
sirkkur2 = r2 + g2 + b2 + s2
sirkkur3 = r3 + g3 + b3 + s3

print "_"*19
print "| Integrert fluks |"
print "|%s|" %("-"*17)
print "| %.11f |\n| %.9f |\n| %.10f |" %(sirkkur1, sirkkur2, sirkkur3)

print "\n"

print "_"*63
print "| Integrert fluks |  Side 1  |  Side 2  |  Side 3  |  Side 4  |"
print "|%s+%s+%s+%s+%s|" %("-"*17, "-"*10, "-"*10, "-"*10, "-"*10)
print "|   Rektangel 1   | %.3f | %.2f | %.2f | %.1f |" %(r1, g1, b1, s1)
print "|   Rektangel 2   | %.2f | %.2f | %.2f | %.1f |" %(r2, g2, b2, s2)
print "|   Rektangel 3   | %.3f | %.3f | %.4f | %.2f |" %(r3, g3, b3, s3)

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig2_g.py
___________________
| Integrert fluks |
|-----------------|
| 104.67598610377 |
| -6366.877538290 |
| -104.3452178478 |


_______________________________________________________________
| Integrert fluks |  Side 1  |  Side 2  |  Side 3  |  Side 4  |
|-----------------+----------+----------+----------+----------|
|   Rektangel 1   | 1556.868 | 19780.11 | -2059.68 | -19172.6 |
|   Rektangel 2   | -5187.56 | 13131.81 | -4074.05 | -10237.1 |
|   Rektangel 3   | -195.570 | 1397.099 | 284.9436 | -1590.82 |

"""