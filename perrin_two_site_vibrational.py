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
bias = 0.15 #in eV
capacitive = 0.0

electron_phonon_coupling = .015*100
phonon_energy = 0.015 # just a low number.

number_of_phonons = 5 
perturbation_expansion_order = 15

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

beta = 250.0

epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);

print "Calculating transmission for order %d" % (perturbation_expansion_order)

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

ret_gf, ad_gf = calculation.singleparticlebackground(0) 



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

minimum = 1.2 * np.min([ np.min(transmission), np.min(old_transmission)])
maximum = 1.2 * np.max([ np.max(transmission), np.max(old_transmission)])

mit = np.min(transmission)
mat = np.max(transmission)
mot = np.min(old_transmission)
maot = np.max(old_transmission)

print "New %.3e < T < %.3e" % (mit, mat)
print "Old %.3e < T < %.3e" % (mot, maot)
if mit > 0 and mat > 0 and mot > 0 and maot > 0:
    plt.semilogy(epsilon, transmission, 'r-', label="vibrational expansion, order=%d" % (perturbation_expansion_order))    
    plt.semilogy(epsilon, old_transmission, 'k--', label="ignore vibrations")   
else:
        plt.plot(epsilon, transmission, 'r-', label="vibrational expansion, order=%d" % (perturbation_expansion_order))    
        plt.plot(epsilon, old_transmission, 'k--', label="ignore vibrations")   
plt.legend()

title = "Transmission, $\\lambda=%.3f" % electron_phonon_coupling
xlabel = "Energy $\\epsilon$"
ylabel = "Transmission"

plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

final_title = "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
    alpha, tau, gamma, levels, bias, beta, capacitive)

print final_title
#plt.title( final_title, fontsize=15)      


non_int = hamiltonian + tunnel

values,_ = np.linalg.eig(non_int)
print values

height = values*0
if plotting_mode == 0 or plotting_mode == 2:
    height += 0.5*maximum

plt.plot(values, height, 'ko') 


if plotting_mode == 2:
    plt.savefig('perrin_two_site_vibrational.svg')
else:    
    plt.show()

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)