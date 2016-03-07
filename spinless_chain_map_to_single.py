#Numeric Python. Should be familiar.
import numpy as np
#Plotting tools
import matplotlib
import matplotlib.pyplot as plt
#Our calculation module
from igf import *
#sys. sys.stderr is a nice text stream.
import sys as sys
#Command line arguments.
import argparse as argparse  
#Time, because I want to see the total calculation time in ms.
import time
#Start the clock.
global_time_start = time.time()
print "Python version is %s.%s.%s., should be >2.7.10 for us. \n" % (sys.version_info[0],sys.version_info[1],sys.version_info[2])
#two level system, to start with
#this is the simplest couloumb-blockade system around

def eigenbasis( matrix, eigenvectors ):
    new_matrix = np.dot( np.dot( np.linalg.inv(eigenvectors), matrix), eigenvectors)
    #new_matrix[ new_matrix < 1e-5] = 0
    
    return new_matrix
# Normally I use argparse for this, but due to sverely limited computational 
#   resources, i.e. a computer that should've been retired years ago, I am
#   working online (Get Data Joy) and can't use commandline arguments.   
plotting_mode = 2
chain_length = 3

capacitive_strength = 0.35
tunnel_strength = -0.5
epsilon_gap = 0.15

epsilon_left = -1
epsilon_right = 3
resolution = 2500

#Inverse temperature (units same as the others)
param_beta = 0.05 * 35

dot_levels = chain_length #All at the same energy; this is a chain.

param_u = np.zeros((dot_levels, dot_levels))
param_tau = np.zeros((dot_levels, dot_levels))

for i in range(0, dot_levels): 
    if i + 1 < chain_length:
        l = i
        r = i+1
        for dl in range(0,1):
            for dr in range(0,1):
                param_tau[l+dl][r+dr] = tunnel_strength;
                param_tau[r+dr][l+dl] = tunnel_strength; 
                param_u[l+dl][r+dr] = capacitive_strength;
                param_u[r+dr][l+dl] = capacitive_strength; 

param_epsilon = np.diag( np.ones((dot_levels)))
#This makes coupling to the leads
#comparable to the coupling between levels. 
gamma_strength = 0.05 

param_gamma_left = np.zeros((dot_levels,dot_levels))
param_gamma_right = np.zeros((dot_levels,dot_levels))

param_gamma_left[0][0] = gamma_strength;
param_gamma_right[epsilon_gap][epsilon_gap] = gamma_strength;


effective_hamiltonian = param_epsilon + param_tau;

print "\nEffective H: \n", effective_hamiltonian

eigenvalues, eigenvectors = np.linalg.eig(effective_hamiltonian)

print effective_hamiltonian
print eigenvalues

new_hamil = eigenbasis( effective_hamiltonian, eigenvectors)
new_epsilon = eigenbasis( param_epsilon, eigenvectors)
new_tau= eigenbasis( param_tau, eigenvectors)
new_u = eigenbasis( param_u, eigenvectors)

new_gamma_left = eigenbasis( param_gamma_left, eigenvectors)
new_gamma_right = eigenbasis( param_gamma_right, eigenvectors)

print "\nEigen H: \n", new_hamil
print "\nEigen Gamma Left: \n", new_gamma_left
print "\nEigen Gamma Right: \n", new_gamma_right



calculation = igfwl(
    param_epsilon, 
    param_tau,
    param_u, 
    param_gamma_left,
    param_gamma_right, 
    param_beta
)
new_calculation = igfwl(
    new_epsilon, 
    new_tau,
    new_u, 
    new_gamma_left,
    new_gamma_right, 
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
    
    transmission = epsilon*0
    
    for i in calculation.generate_superset(0):
        transmission += calculation.transport_channel(i, epsilon)
    
    plt.plot(epsilon, transmission, 'g-', label="Original") 
    
    transmission = epsilon*0
    
    for i in new_calculation.generate_superset(0):
        transmission += new_calculation.transport_channel(i, epsilon)
    
    plt.plot(epsilon, transmission, 'r--', label="Mapped (NI)") 
     
    maximum = 0.02
    if np.max(transmission) > maximum:
        maximum = np.max(transmission)*1.25
    
    plt.ylim([0, maximum])
    
    title = "Chain Transmission"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Transmission"
    
elif plotting_mode == 1 or plotting_mode == 3:
    
    
    spectral = epsilon*0 
    
    for i in calculation.generate_superset(0):
        spectral += calculation.spectral_channel(i, epsilon)
    
    plt.plot(epsilon, spectral, 'g-', label="Original")  
    
    spectral = epsilon*0 
    
    for i in calculation.generate_superset(0):
        spectral += calculation.spectral_channel(i, epsilon)
    
    plt.plot(epsilon, spectral, 'g-', label="Mapped (NI)")  
    
    
    title = "Chain Spectral"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Spectral"

plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "%d-%s: $\\beta=%.3f$, $\\epsilon_0=%.3f$, $\\Gamma=%.3f$, $\\tau=%.3f$ , $U=%.3f$" % (chain_length, title, calculation.beta,
    epsilon_gap, gamma_strength, tunnel_strength, capacitive_strength), fontsize=15)     
plt.legend()


if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('spinless_chain_map.svg')
else:    
    plt.show()

global_time_end = time.time ()

print "\n Time spent %.6f seconds in mode %d. \n " % (global_time_end - global_time_start, plotting_mode)