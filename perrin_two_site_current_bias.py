import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
global_time_start = time.time()

plotting_mode = 0

alpha = 0.74
tau = 0.0241
gamma = 0.0102
levels = -0.0 #Fig 4b, Perrin(2014)
capacitive = 0.45 

tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

# left, right are now +- eV/2, see Fig 4b in Perrin(2014)
epsilon_res = 100

bias_left = -1.0
bias_right = 1.0
bias_res = 100

interaction = np.zeros((2,2))
interaction[0][1] = capacitive
interaction[1][0] = capacitive


gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

beta = 0.05 * 150

biaswindow = np.linspace(bias_left, bias_right, bias_res)
    
current = []

for bias in biaswindow:
    hamiltonian = np.zeros((2,2))
    
    hamiltonian[0][0] = levels + 0.5 * alpha * bias
    hamiltonian[1][1] = levels - 0.5 * alpha * bias
    
    calculation = igfwl(
        hamiltonian, 
        tunnel,
        interaction, 
        gamma_left,
        gamma_right, 
        beta
    )
    
    epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
    
    #It is unfeasible to plot all the channels. Sum them up!
    
    transmission = epsilon*0
    
    for i in calculation.generate_superset(0):
        transmission += calculation.transport_channel(i, epsilon)
    
    current.append( [ np.trapz(transmission, epsilon) ])


current = np.array(current)
minimum = 1.2 * np.min(current)
maximum = 1.2 * np.max(current)



plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')
 
plt.plot(biaswindow, current, 'g-')   
 

title = "Current versus bias"
xlabel = "Bias $V_b$"
ylabel = "$I(V_b)$"
 
plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
    alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=15)     
#plt.legend()


if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('perrin_two_site.svg')
else:    
    plt.show()

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)