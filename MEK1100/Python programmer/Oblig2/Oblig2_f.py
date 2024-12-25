#coding=utf-8
from Oblig2_e import *
from numpy import sum

r1 = sum(u[159, 34:70])*0.5; g1 = sum(v[159:169, 69])*0.5; b1 = - sum(u[169, 34:70])*0.5; s1 = - sum(v[159:169, 34])*0.5
r2 = sum(u[ 84, 34:70])*0.5; g2 = sum(v[ 84: 99, 69])*0.5; b2 = - sum(u[ 99, 34:70])*0.5; s2 = - sum(v[ 84: 99, 34])*0.5
r3 = sum(u[ 49, 34:70])*0.5; g3 = sum(v[ 49: 59, 69])*0.5; b3 = - sum(u[ 59, 34:70])*0.5; s3 = - sum(v[ 49: 59, 34])*0.5

sirkkur1 = r1 + g1 + b1 + s1
sirkkur2 = r2 + g2 + b2 + s2
sirkkur3 = r3 + g3 + b3 + s3

sirkfla1 = sum(curlz[159:170, 34:70])*0.5
sirkfla2 = sum(curlz[84:100, 34:70])*0.5
sirkfla3 = sum(curlz[49:60, 34:70])*0.5

print "_"*63
print "| Sirkulasjon |  rektangel 1  |  rektangel 2  |  rektangel 3  |"
print "|%s+%s+%s+%s|" %("-"*13, "-"*15, "-"*15, "-"*15)
print "|Kurveintegral|  %.6f  |  %.4f  |  %.7f  |" %(sirkkur1, sirkkur2, sirkkur3)
print "|Flateintegral|  %.6f  |  %.4f  |  %.7f  |" %(sirkfla1, sirkfla2, sirkfla3)
print "\n"

print "_"*59
print "| Sirkulasjon |  Side 1  |  Side 2  |  Side 3  |  Side 4  |"
print "|%s+%s+%s+%s+%s|" %("-"*13, "-"*10, "-"*10, "-"*10, "-"*10)
print "| Rektangel 1 | %.2f | %.4f | %.1f | %.4f |" %(r1, g1, b1, s1)
print "| Rektangel 2 | %.4f | %.4f | %.1f | %.3f |" %(r2, g2, b2, s2)
print "| Rektangel 3 | %.3f | %.4f | %.2f | %.5f |" %(r3, g3, b3, s3)

"""
C:\Users\Torstein\Documents\UiO\Mek1100\Python programmer>Oblig2_f.py
| Sirkulasjon |  rektangel 1  |  rektangel 2  |  rektangel 3  |
|-------------+---------------+---------------+---------------|
|Kurveintegral|  2635.667121  |  -60819.0999  |  -21.0157324  |
|Flateintegral|  2621.558696  |  -61482.5410  |  -12.2143339  |


___________________________________________________________
| Sirkulasjon |  Side 1  |  Side 2  |  Side 3  |  Side 4  |
|-------------+----------+----------+----------+----------|
| Rektangel 1 | 70100.52 | 235.1629 | -68332.9 | 632.8364 |
| Rektangel 2 | 198.4756 | 462.2918 | -61243.5 | -236.403 |
| Rektangel 3 | 5133.348 | 185.4207 | -5410.04 | 70.25540 |

"""


