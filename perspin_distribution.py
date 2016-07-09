import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
#Command line arguments.
from scipy.optimize import minimize
import argparse as argparse  
import time
global_time_start = time.time()


parser  = argparse.ArgumentParser(prog="perspin transport map",
  description = "Compare spinfull/spinless transmissions.")  
  
parser.add_argument(
    '-k',
    '--ksi',
    help='Onsite interaction',
    action='store',
    type = float,
    default = 1.0
)   
parser.add_argument(
    '-u',
    '--capacitive',
    help='Onsite interaction',
    action='store',
    type = float,
    default = 0.5
)   


#change these numbers based on visual inspection
epsilon_left = -.5
epsilon_right = .5
epsilon_res = 250
epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res); 

alpha = 0.5
tau = 0.02
gamma = 0.01
bias = 0.25

args    = parser.parse_args()  
capacitive = args.capacitive
ksi = args.ksi

levels = np.linspace(-0.5, -1e-6, 25)

E_zero_zero_zero_zero   = 0 * levels 
E_zero_zero_zero_one    = levels - 0.5 * alpha * bias
E_zero_zero_one_zero    = levels - 0.5 * alpha * bias
E_zero_zero_one_one     = 2*(levels - 0.5 * alpha * bias) + ksi * capacitive

E_zero_one_zero_zero    = levels + 0.5 * alpha * bias
E_zero_one_zero_one     = 2 * levels + capacitive
E_zero_one_one_zero     = 2* levels + capacitive
E_zero_one_one_one      = 3 * levels - 0.5 * alpha * bias + capacitive* 2 + ksi * capacitive

E_one_zero_zero_zero    = levels + 0.5 * alpha * bias
E_one_zero_zero_one     = 2 * levels + capacitive
E_one_zero_one_zero     = 2 * levels  + capacitive
E_one_zero_one_one      = 3 * levels - 0.5 * alpha * bias + capacitive* 2 + ksi * capacitive

E_one_one_zero_zero     = 2 * levels + 0.5 * alpha * bias * 2 + ksi * capacitive
E_one_one_zero_one      = 3 * levels + 0.5 * alpha * bias + capacitive * 2+ ksi*capacitive
E_one_one_one_zero      = 3 * levels + 0.5 * alpha * bias + capacitive * 2+ ksi*capacitive
E_one_one_one_one       = 4 * levels + 2 * capacitive + 2 * capacitive * ksi + 2 * capacitive
    
plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
 
 
title = " $\\alpha=%.3f$, $V=%.3f$, $U=%.3f$, $\\xi=%.3f$" % (alpha, bias, capacitive, ksi)
xlabel = "Zero-bias level $\\epsilon_0$"
ylabel = "Energy [eV]"
plt.rc('font', family='serif')

plt.plot(levels, E_zero_zero_zero_zero, 'ro', label='$E_{0000}$')     
plt.plot(levels, E_zero_zero_zero_one, 'rd', label='$E_{0001}$')     
plt.plot(levels, E_zero_zero_one_zero, 'r^', label='$E_{0010}$') 
plt.plot(levels, E_zero_zero_one_one, 'rv', label='$E_{0011}$')         

plt.plot(levels, E_zero_one_zero_zero, 'go', label='$E_{0100}$')     
plt.plot(levels, E_zero_one_zero_one, 'gd', label='$E_{0101}$')     
plt.plot(levels, E_zero_one_one_zero, 'g^', label='$E_{0110}$') 
plt.plot(levels, E_zero_one_one_one, 'gv', label='$E_{0111}$')         

plt.plot(levels, E_one_zero_zero_zero, 'bo', label='$E_{1000}$')     
plt.plot(levels, E_one_zero_zero_one, 'bd', label='$E_{1001}$')     
plt.plot(levels, E_one_zero_one_zero, 'b^', label='$E_{1010}$') 
plt.plot(levels, E_one_zero_one_one, 'bv', label='$E_{1011}$')         

plt.plot(levels, E_one_one_zero_zero, 'mo', label='$E_{1100}$')     
plt.plot(levels, E_one_one_zero_one, 'md', label='$E_{1101}$')     
plt.plot(levels, E_one_one_one_zero, 'm^', label='$E_{1110}$') 
plt.plot(levels, E_one_one_one_one, 'mv', label='$E_{1111}$')         


plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

plt.title( "%s" % (title), fontsize=15)     
plt.legend(loc='upper left')

plt.savefig('perspin_distribution.pdf') 