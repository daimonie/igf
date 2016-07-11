import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
#Command line arguments.
from scipy.optimize import minimize
import argparse as argparse  
from scipy.optimize import curve_fit 



capacitive = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, .35, .40, .45, .50])

#parameters
#alpha   = 0.50
#tau     = 0.02
#gamma   = 0.01
#zeta    = 1.00
#ksi     = 1.00
#NB: Current is in nA

#spinless level -0
current_spinless_zero = np.array([6.382e2, 1.19e2, 4.872e1, 2.579e1, 1.603e1, 1.203e01, 9.19, 8.643, 7.801, 7.227, 6.802])
#spinless level -u
current_spinless_cap = np.array([6.382e2, 1.238e2, 5.45e1, 3.518e1, 2.809e1, 2.46e1, 2.16e1, 2.135e1, 2.049e1, 1.986e1, 1.941e1])
 


##plotting
plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
 
 
#title = " $\\alpha=%.3f$, $V=%.3f$, $U=%.3f$" % (alpha, bias, capacitive)
xlabel = "Capacitive Interaction strength $U$"
ylabel = "Maximum current $I$ [nA]"
plt.rc('font', family='serif')

plt.semilogy(capacitive[0:len(current_spinless_zero)], current_spinless_zero, 'rd', label='Spinless $\\epsilon_0=0$') 

polyfitter = np.polyfit( capacitive, current_spinless_zero, 2)

testfunc = lambda u, a, b, c, d, e: a * np.exp(-c*u) + b * np.exp(-d*u) + e

popt, pcov =  curve_fit(testfunc, capacitive, current_spinless_zero)
print "spinless", popt
plt.semilogy(capacitive, testfunc(capacitive, *popt), 'r--', label='Spinless fit $ae^{-bu}+ce^{-du}+e$') 

popt, pcov =  curve_fit(testfunc, capacitive, current_spinless_cap)
print "spinless", popt
plt.semilogy(capacitive, testfunc(capacitive, *popt), 'g--', label='Spinfull fit $ae^{-bu}+ce^{-du}+e$') 


plt.semilogy(capacitive[0:len(current_spinless_cap)], current_spinless_cap, 'g^', label='Spinless $\\epsilon_0=-U$')       

plt.xticks(capacitive)

plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

#plt.title( "%s" % (title), fontsize=15)     
plt.legend(loc='upper right')

plt.savefig('imax.pdf') 