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
levels = -0.25
bias = 0.00 #in eV
capacitive = 0.2

hamiltonian = np.zeros((2,2))

hamiltonian[0][0] = levels + 0.5 * alpha * bias
hamiltonian[1][1] = levels - 0.5 * alpha * bias

tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

#change these numbers based on visual inspection
epsilon_left = -0.4
epsilon_right = -0.15
epsilon_res = 10000

interaction = np.zeros((2,2))
interaction[0][1] = capacitive
interaction[1][0] = capacitive


gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

beta = 0.05 * 50*5


calculation = igfwl(
    hamiltonian, 
    tunnel,
    interaction, 
    gamma_left,
    gamma_right, 
    beta
)


epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);


plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

maximum = 0.2
if plotting_mode == 0 or plotting_mode == 2:
    
    #It is unfeasible to plot all the channels. Sum them up!
    
    transmission = epsilon*0
    
    chances = calculation.distribution()
    
    for k in calculation.generate_superset(0):
        this_channel = calculation.transport_channel(k, epsilon)
        print "P[%d]=%.3f,\t%.3f\t%.3f" % (k, chances[k], np.min(this_channel), np.max(this_channel))
        transmission += this_channel
    print "Total\t%2.3f\t%2.3f" % (np.min(transmission), np.max(transmission))
    
    transmission /= .355
    maximum = 1.2 * np.max(transmission)
    #plt.semilogy(epsilon, transmission, 'g-')   
    title = "Transmission"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Transmission"
    
    plt.plot(epsilon, transmission, 'g-',label="Many-body %s" % title)   
     
    
    
elif plotting_mode == 1 or plotting_mode == 3:
    
    
    spectral = epsilon*0 
    
    for k in calculation.generate_superset(0):
        spectral += calculation.spectral_channel(k, epsilon)
    
    maximum = 1.2 * np.max(spectral)
    plt.plot(epsilon, spectral, 'g-', label="sum")  
    
    title = "Spectral"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Spectral"

plt.ylim([0, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
    alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=15)     
plt.legend()


non_int = hamiltonian + tunnel

values,_ = np.linalg.eig(non_int)
print values

height = values*0
if plotting_mode == 0 or plotting_mode == 2:
    height += np.average(transmission)
elif plotting_mode == 1 or plotting_mode == 3:
    height += np.average(spectral)

plt.plot(values, height, 'ko') 


if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('perrin_two_site.svg')
else:    
    plt.show()

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
