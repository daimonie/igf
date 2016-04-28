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
tau = 0.005
gamma = 0.01
bias = 0.25*1
capacitive = .100*1
ksi = 2.0


#levels set to zero-
levels = -0.05
levels= -capacitive

hamiltonian = np.zeros((4,4))

hamiltonian[0][0] = levels + 0.5 * alpha * bias
hamiltonian[1][1] = levels + 0.5 * alpha * bias
hamiltonian[2][2] = levels - 0.5 * alpha * bias
hamiltonian[3][3] = levels - 0.5 * alpha * bias

tunnel = np.zeros((4,4))
tunnel[0][2] = tau 
tunnel[0][3] = tau 

tunnel[1][2] = tau 
tunnel[1][3] = tau 

tunnel[2][0] = tau 
tunnel[2][1] = tau

tunnel[3][0] = tau 
tunnel[3][1] = tau 
 
#change these numbers based on visual inspection
epsilon_left = -.6
epsilon_right = .4
epsilon_res = 10000

interaction = np.zeros((4,4))

interaction[0][1] = ksi * capacitive
interaction[1][0] = ksi * capacitive

interaction[2][3] = ksi * capacitive
interaction[3][2] = ksi * capacitive

interaction[0][3] = capacitive
interaction[3][0] = capacitive

interaction[0][2] = capacitive
interaction[2][0] = capacitive

interaction[1][3] = capacitive
interaction[3][1] = capacitive

interaction[1][2] = capacitive
interaction[2][1] = capacitive


##print interaction
#sys.exit(0)
gamma_left = np.zeros((4,4))
gamma_left[0][0] = gamma
gamma_left[1][1] = gamma

gamma_right = np.zeros((4,4))
gamma_right[2][2] = gamma
gamma_right[3][3] = gamma

beta = 250.0*1


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


fig = plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
fig.subplots_adjust(left=0.17)

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
     
    maximum = 1.2 * np.max(transmission*scaler)
    title = "Transmission"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Transmission"
    
    minimum = 1e-6
    #plt.semilogy(epsilon, transmission , 'r--',label="scaled(Many-body %s)" % title)   
    #plt.semilogy(epsilon, transmission * scaler, 'g-',label="Many-body %s" % title)   
    #plt.plot(epsilon, transmission , 'r--',label="scaled(Many-body %s)" % title)   
    plt.plot(epsilon, transmission * scaler, 'g-',label="Many-body %s" % title)   
     
    
    
elif plotting_mode == 1 or plotting_mode == 3:
    
    
    spectral = epsilon*0 
    
    for k in calculation.generate_superset(0):
        spectral += calculation.spectral_channel(k, epsilon)
    
    maximum = 1.2 * np.max(spectral)
    plt.plot(epsilon, spectral, 'g-', label="sum")  
    
    title = "Spectral"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Spectral"

if maximum < minimum:
    maximum = 1.0
    
plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "Per.Spin [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
    alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=15)     
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


if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('perspin_two_site.svg')
else:    
    plt.savefig('perspin_two_site.png')

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
