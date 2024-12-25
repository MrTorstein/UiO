from Oblig1_e import A, Euler #Importerer funksjonene A og Euler fra forrige fil

t, s, v, a = Euler(A)

print "Loperen bruker %.2fsek paa 100m" %t[-1]	#Henter ut siste ledd i listen t, som er der loperen naar 100m merket


"""
C:\Users\Torstein\Documents\UiO\Fys-Mek1110\Python Programmer>Oblig1_f.py
Loperen bruker 6.79sek paa 100m

"""