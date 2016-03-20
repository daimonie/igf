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

number_of_phonons = 5 
perturbation_expansion_order = 5

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
    electron_phonon_coupling,
    number_of_phonons,
    perturbation_expansion_order
) 
epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);

transmission = calculation.full_transport(epsilon)




old_calculation = igfwl(
    hamiltonian, 
    tunnel,
    interaction, 
    gamma_left,
    gamma_right, 
    beta
)


old_transmission = epsilon*0

for i in old_calculation.generate_superset(0):
    old_transmission += old_calculation.transport_channel(i, epsilon)

plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

maximum = 1.0

minimum = 1.2 * np.min(transmission)
maximum = 1.2 * np.max(transmission)

plt.plot(epsilon, transmission, 'g-')   
plt.plot(epsilon, old_transmission, 'r--')       
    

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