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
bias = 0.25*2
capacitive = .15

#levels set to zero-
#levels = -0.05
levels=-capacitive


gamma *= 3

hamiltonian = np.zeros((2,2))

hamiltonian[0][0] = levels + 0.5 * alpha * bias
hamiltonian[1][1] = levels - 0.5 * alpha * bias

tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

#change these numbers based on visual inspection
epsilon_left = -.6
epsilon_right = .4
epsilon_res = 10000

interaction = np.zeros((2,2))
interaction[0][1] = capacitive
interaction[1][0] = capacitive


gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

beta = 250.0*1



calculation = igfwl(
    hamiltonian, 
    tunnel,
    interaction, 
    gamma_left,
    gamma_right, 
    beta
)

print "Chance overview:"
p = calculation.distribution()
for i in calculation.generate_superset(0):
    print "%d\t%.3e" % (i, p[i])


epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);


plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

minimum = 0.0
maximum = 0.2
if plotting_mode == 0 or plotting_mode == 2:
    
    #It is unfeasible to plot all the channels. Sum them up!
    
    transmission = calculation.full_transmission(epsilon)
    
    scaler = calculation.scaler()
        
    print "Total\t%2.3f\t%2.3f" % (np.min(transmission), np.max(transmission))
    print "Scaler:\t%2.3f" % scaler
     
    maximum = 1.2 * np.max(transmission)
    title = "Transmission $T(\\epsilon)$"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Transmission"
    
    #plt.semilogy(epsilon, transmission , 'r--',label="scaled(Many-body %s)" % title)   
    #plt.semilogy(epsilon, transmission * scaler, 'g-',label="Many-body %s" % title)   
    #plt.plot(epsilon, transmission , 'r--',label="scaled(Many-body %s)" % title)   
    plt.plot(epsilon, transmission * scaler, 'g-', linewidth=5, label="Interacting transmission")   
    
    if capacitive > 0.00:
        non_interacting_calculation = igfwl(
            hamiltonian, 
            tunnel,
            interaction*0, 
            gamma_left,
            gamma_right, 
            beta
        )
        
        non_interacting_transmission = non_interacting_calculation.full_transmission(epsilon)  
        plt.plot(epsilon, non_interacting_transmission, 'r-', linewidth=5 ,label="Non-Interacting")   
        plt.legend()
        
        minimum = 0.00
        maximum = non_interacting_transmission.max()
        
    
elif plotting_mode == 1 or plotting_mode == 3:
    
    
    spectral = epsilon*0 
    
    for k in calculation.generate_superset(0):
        spectral += calculation.spectral_channel(k, epsilon)
    
    maximum = 1.2 * np.max(spectral)
    plt.plot(epsilon, spectral, 'g-', label="sum")  
    
    title = "Spectral"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Spectral"

plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
    alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=15)     
#plt.legend()


non_int = hamiltonian + tunnel

values,_ = np.linalg.eig(non_int)
print "eigs:", values

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
