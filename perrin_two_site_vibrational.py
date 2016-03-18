import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
global_time_start = time.time()

plotting_mode = 2

alpha = 0.74
tau = 0.0241
gamma = 0.0102
levels = -0.25
bias = 1.0 #in eV
capacitive = 0.0

electron_phonon_coupling = 1.5
phonon_energy = 0.015 # just a low number.

global_max_phonons = 5*15 #if set too high, double factorials leads to Inf.

cutoff_chance = 1e-4

hamiltonian = np.zeros((2,2))

hamiltonian[0][0] = levels + 0.5 * alpha * bias
hamiltonian[1][1] = levels - 0.5 * alpha * bias

tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

epsilon_left = - 2.0
epsilon_right = 2.0
epsilon_res = 1000

interaction = np.zeros((2,2))
interaction[0][1] = capacitive
interaction[1][0] = capacitive


gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

beta = 0.05 * 150


calculation = igfwl_vibrational(
    hamiltonian, 
    tunnel,
    interaction, 
    gamma_left,
    gamma_right, 
    beta,
    phonon_energy,
    electron_phonon_coupling    
)
calculation.cutoff_chance = cutoff_chance

max_phonon_number = 0

chances = []

for i in range(global_max_phonons):
    chance = calculation.chance_phonons(i, i)
    chances.append(chance)
    
    if chance < calculation.cutoff_chance:
        max_phonon_number = i
        break

summation_factors = np.zeros((max_phonon_number,max_phonon_number,max_phonon_number,max_phonon_number))    
print "Maximum phonon number by boltzman is %d" % max_phonon_number

if max_phonon_number > global_max_phonons:
    max_phonon_number = global_max_phonons
    print "Maximum phonon number exceeds calculational options; setting to %d." % max_phonon_number

significant_terms = []    
for n in range(max_phonon_number):
    for m in range(max_phonon_number):
        if chances[m] > calculation.cutoff_chance:
            for x in range(max_phonon_number): 
                if (calculation.pe * calculation.pec)**x > calculation.cutoff_chance:
                    for y in range(x):  
                        prefactor = calculation.sum_factor(n,m,x,y) * chances[m]
                        
                        if prefactor < calculation.cutoff_chance:
                            break
                        elif prefactor > calculation.cutoff_chance:
                            #print "p[%d,%d,%d,%d] = %.3e" % (n,m,x,y, prefactor)   
                            significant_terms.append( (n,m,x,y, prefactor))
                    
print "Found %d significant terms (under max. phonon constraint)" % len(significant_terms)

epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);


plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

maximum = 0.2

transmission = epsilon*0

for i in calculation.generate_superset(0):
    transmission += calculation.transport_channel_vibrational(i, epsilon, significant_terms)

minimum = 1.2 * np.min(transmission)
maximum = 1.2 * np.max(transmission)
#plt.plot(epsilon, transmission, 'g-')   
plt.semilogy(epsilon, transmission, 'g-')   
    

title = "Transmission, $\\lambda=%.3f" % electron_phonon_coupling
xlabel = "Energy $\\epsilon$"
ylabel = "Transmission"

plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

final_title = "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$,max=%.3e" % (title,
    alpha, tau, gamma, levels, bias, beta, capacitive, np.max(transmission))

print final_title
plt.title( final_title, fontsize=15)     
#plt.legend()


non_int = hamiltonian + tunnel

values,_ = np.linalg.eig(non_int)
print values

height = values*0
if plotting_mode == 0 or plotting_mode == 2:
    height += np.average(transmission) 

plt.plot(values, height, 'ko') 


if plotting_mode == 2:
    plt.savefig('perrin_two_site_vibrational.svg')
else:    
    plt.show()

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)