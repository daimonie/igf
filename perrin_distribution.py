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

E_zero_zero = 0 * levels
E_one_zero  = levels + 0.5 * alpha * bias
E_zero_one  = levels - 0.5 * alpha * bias
E_one_one   = 2 * levels + capacitive
    
plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
 
 
title = " $\\alpha=%.3f$, $V=%.3f$, $U=%.3f$" % (alpha, bias, capacitive)
xlabel = "Zero-bias level $\\epsilon_0$"
ylabel = "Energy [eV]"
plt.rc('font', family='serif')

plt.plot(levels, E_zero_zero, 'ro', label='$E_{00}$')     
plt.plot(levels, E_one_zero, 'gd', label='$E_{10}$')     
plt.plot(levels, E_zero_one, 'b^', label='$E_{01}$') 
plt.plot(levels, E_one_one, 'mv', label='$E_{11}$')         


plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

plt.title( "%s" % (title), fontsize=15)     
plt.legend(loc='upper left')

plt.savefig('perrin_distribution.pdf') 