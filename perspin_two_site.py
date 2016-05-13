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

beta        = 250.00  
gamma       = 0.010
tau         = 0.020
alpha       = 0.300
capacitive  = 0.117
alpha       = 0.25
levels      = -0.122
bias        = 0.00

ksi = 2.0*1e5
zeta = 1.0

hamiltonian = np.zeros((4,4))

hamiltonian[0][0] = levels + 0.5 * alpha * bias
hamiltonian[1][1] = levels + 0.5 * alpha * bias
hamiltonian[2][2] = levels - 0.5 * alpha * bias
hamiltonian[3][3] = levels - 0.5 * alpha * bias

tunnel = np.zeros((4,4))
tunnel[0][2] = tau  
tunnel[1][3] = tau  

#right-left
tunnel[2][0] = tau  
tunnel[3][1] = tau  

#change these numbers based on visual inspection
epsilon_left = -.3
epsilon_right = .2
epsilon_res = 10000

interaction = np.zeros((4,4))

interaction[0][1] = ksi * capacitive
interaction[1][0] = ksi * capacitive

interaction[2][3] = ksi * capacitive
interaction[3][2] = ksi * capacitive

interaction[0][3] = zeta*capacitive
interaction[3][0] = zeta*capacitive

interaction[0][2] = zeta*capacitive
interaction[2][0] = zeta*capacitive

interaction[1][3] = zeta*capacitive
interaction[3][1] = zeta*capacitive

interaction[1][2] = zeta*capacitive
interaction[2][1] = zeta*capacitive

##print interaction
#sys.exit(0)
gamma_left = np.zeros((4,4))
gamma_left[0][0] = gamma
gamma_left[1][1] = gamma

gamma_right = np.zeros((4,4))
gamma_right[2][2] = gamma
gamma_right[3][3] = gamma
 

calculation = igfwl(
    hamiltonian, 
    tunnel,
    interaction, 
    gamma_left,
    gamma_right, 
    beta
)
print hamiltonian
print interaction
print "Chance overview:"
p = calculation.distribution()
for i in calculation.generate_superset(0):
    x = calculation.ket(i)
    print "%d\t%.3e\t%d\t%d\t%d\t%d" % (i, p[i], x[0], x[1], x[2], x[3])
epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);


#fig = plt.figure(figsize=(10, 10), dpi=1080)
fig = plt.figure()
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
fig.subplots_adjust(left=0.17)
 
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

minimum = 0.0
maximum = 2.0 

transmission = calculation.full_transmission(epsilon)  
  
spinless_hamiltonian = np.zeros((2,2))

spinless_hamiltonian[0][0] = levels + 0.5 * alpha * bias
spinless_hamiltonian[1][1] = levels - 0.5 * alpha * bias

spinless_tunnel = np.zeros((2,2))
spinless_tunnel[0][1] = -tau
spinless_tunnel[1][0] = -tau

spinless_interaction = np.zeros((2,2))
spinless_interaction[0][1] = capacitive
spinless_interaction[1][0] = capacitive


spinless_gamma_left = np.zeros((2,2))
spinless_gamma_left[0][0] = gamma

spinless_gamma_right = np.zeros((2,2))
spinless_gamma_right[1][1] = gamma

spinless_calculation = igfwl(
    spinless_hamiltonian, 
    spinless_tunnel,
    spinless_interaction, 
    spinless_gamma_left,
    spinless_gamma_right, 
    beta
)

spinless_transmission = spinless_calculation.full_transmission(epsilon)

maximum = 1.2 * np.max([np.max(transmission),np.max(spinless_transmission)])
   
xlabel = "Energy $\\epsilon$"
ylabel = "Transmission"



minimum = 1e-6
plt.plot(epsilon, transmission, 'g-',label="Many-body (spinful)")   
plt.plot(epsilon, spinless_transmission, 'r--',label="Many-body (spinless)")   
plt.legend()

if maximum < minimum:
    #maximum = 1.0
    minimum = 0.00
plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "$\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=11)     
#plt.legend()


non_int = hamiltonian + tunnel
values,_ = np.linalg.eig(non_int)
print "eigs:", values
height = values*0
if plotting_mode == 0 or plotting_mode == 2:
    height += 0.2 * maximum
elif plotting_mode == 1 or plotting_mode == 3:
    height += 0.2 * maximum

plt.plot(values, height, 'ko') 


spinless_non_int = spinless_hamiltonian + spinless_tunnel
spinless_values,_ = np.linalg.eig(non_int)
print "eigs:", spinless_values
spinless_height = spinless_values*0
if plotting_mode == 0 or plotting_mode == 2:
    spinless_height += 0.2 * maximum
elif plotting_mode == 1 or plotting_mode == 3:
    spinless_height += 0.2 * maximum

plt.plot(spinless_values, spinless_height, 'mo') 


if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('perspin_two_site.svg')
else:    
    plt.show()
    #plt.savefig('perspin_two_site.png')

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
