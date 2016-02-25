import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
import time
global_time_start = time.time()
#two level system, to start with
#this is the simplest couloumb-blockade system around

#units are arbitrary. 
print "Python version is %s.%s.%s., should be >2.7.10 for us. \n" % (sys.version_info[0],sys.version_info[1],sys.version_info[2])


# Normally I use argparse for this, but due to sverely limited computational 
#   resources, i.e. a computer that should've been retired years ago, I a
#   working online (Get Data Joy) and can't use commandline arguments.   
plotting_mode = 2
chain_length = 4

dot_levels = 2* chain_length #All at the same energy; this is a chain.

capacitive_strength = 00.
tunnel_strength = 0.5
epsilon_gap = .15

param_u = np.zeros((dot_levels, dot_levels))
param_tau = np.zeros((dot_levels, dot_levels))

for i in range(0, dot_levels):
    if i%2==0:
        l = i
        r = l + 1
        
        param_u[l][r] = capacitive_strength;
        param_u[r][l] = capacitive_strength;
        
    if i%2==0 and i/2 < chain_length-1:
        l = i
        r = i+2
        for dl in range(0,2):
            for dr in range(0,2):
                param_tau[l+dl][r+dr] = tunnel_strength;
                param_tau[r+dr][l+dl] = tunnel_strength;

print param_tau
print param_u
param_epsilon = np.diag( np.ones((dot_levels)))
epsilon_left = -5.0
epsilon_right = 5.0
resolution = 100

#This makes coupling to the leads
#comparable to the coupling between levels. 
gamma_strength = 0.05 

param_gamma_left = np.zeros((dot_levels,dot_levels))
param_gamma_right = np.zeros((dot_levels,dot_levels))

param_gamma_left[0][0] = gamma_strength;
param_gamma_right[epsilon_gap][epsilon_gap] = gamma_strength;

#Inverse temperature (units same as the others)
param_beta = 0.05 * 50

calculation = igfwl(
    param_epsilon, 
    param_tau,
    param_u, 
    param_gamma_left,
    param_gamma_right, 
    param_beta
)

epsilon = np.linspace(epsilon_left, epsilon_right, resolution)
plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

if plotting_mode == 0 or plotting_mode == 2:
    
    #It is unfeasible to plot all the channels. Sum them up!
    
    transmission = calculation.transport_channel(0, epsilon)
    
    for i in calculation.generate_superset(0):
        transmission += calculation.transport_channel(i, epsilon)
    
    plt.plot(epsilon, transmission, 'g-', label="sum")  
    plt.xlabel("energy")
    plt.ylabel("Transmission")
    
    maximum = 0.02
    if np.max(transmission) > maximum:
        maximum = np.max(transmission)*1.25
    
    plt.ylim([0, maximum])
    
    title = "Chain Transmission"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Transmission"
    
elif plotting_mode == 1 or plotting_mode == 3:
    
    
    spectral = calculation.spectral_channel(0, epsilon)
    
    for i in calculation.generate_superset(0):
        spectral += calculation.spectral_channel(i, epsilon)
    
    plt.plot(epsilon, spectral, 'g-', label="sum")  
    
    
    title = "Chain Spectral"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Spectral"

plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "%d-%s: $\\beta=%.3f$, $\\epsilon_0=%.3f$, $\\Gamma=%.3f$, $\\tau=%.3f$ , $U=%.3f$" % (chain_length, title, calculation.beta,
    epsilon_gap, gamma_strength, tunnel_strength, capacitive_strength), fontsize=15)     
plt.legend()

if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('chain.svg')
    
plt.show()

global_time_end = time.time ()

print "\n Time spent %.6f seconds in mode %d. \n " % (global_time_end - global_time_start, plotting_mode)