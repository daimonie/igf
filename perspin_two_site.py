import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator 
from matplotlib.ticker import FormatStrFormatter
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
bias = 0.25*2
capacitive = .15
beta = 250.00

#change these numbers based on visual inspection
epsilon_left = -0.50
epsilon_right = 1.00
epsilon_res = 10000

#levels set to zero-
#levels = -0.05
levels=-capacitive*0


#gamma *= 3

ksi = 3.5
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

print "Chance overview:"
p = calculation.distribution()
for i in calculation.generate_superset(0):
    print "%d\t%.3e" % (i, p[i])


epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);

 
fig, ax = plt.subplots(1, 1, figsize=(25, 15), dpi=1080)
#fig = plt.figure()
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
fig.subplots_adjust(left=0.17)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')


decorative = False
decorative = True

if plotting_mode == 0 or plotting_mode == 2:
    
    #It is unfeasible to plot all the channels. Sum them up!
    
    transmission = calculation.full_transmission(epsilon)
    
    scaler = calculation.scaler()
        
    print "Total\t%2.3f\t%2.3f" % (np.min(transmission), np.max(transmission))
    print "Scaler:\t%2.3f" % scaler
     
    maximum = 1.2 * np.max(transmission)
    title = "T"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Transmission $T(\\epsilon)$"
    
    #plt.semilogy(epsilon, transmission , 'r--',label="scaled(Many-body %s)" % title)   
    #plt.semilogy(epsilon, transmission * scaler, 'g-',label="Many-body %s" % title)   
    #plt.plot(epsilon, transmission , 'r--',label="scaled(Many-body %s)" % title)   
    #plt.plot(epsilon, transmission * scaler, 'g-', linewidth=2, label="Transmission $T(\\epsilon)$")   
    
    if decorative:
        plt.fill(epsilon, transmission * scaler, 'g-', linewidth=2, alpha=0.9, label="Transmission $T(\\epsilon)$")   
    else:
        plt.plot(epsilon, transmission * scaler, 'g-', linewidth=2, label="Transmission $T(\\epsilon)$")   
    minimum = 0.00
    maximum = transmission.max()
        
    
elif plotting_mode == 1 or plotting_mode == 3:
    
    
    spectral = epsilon*0 
    
    for k in calculation.generate_superset(0):
        spectral += calculation.spectral_channel(k, epsilon)
    
    maximum = 1.2 * np.max(spectral)
    plt.plot(epsilon, spectral, 'g-', linewidth=2, label="Spectral $A\\left(\\epsilon\\right)$")  
    
    title = "S"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Spectral"
plt.xlim([-0.50, 1.00])
plt.ylim([-0.001/4.0, 0.27])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30, labelpad=15)

 
plt.xticks(np.array(range(7))*0.25 - 0.50) 
plt.yticks(np.array(range(7))*0.27/6.00)
 
minorLocator1 = AutoMinorLocator(5)
minorLocator2 = AutoMinorLocator(5)

ax.xaxis.set_minor_locator(minorLocator1) 
ax.yaxis.set_minor_locator(minorLocator2) 

ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=20)
plt.tick_params(which='minor', length=10)



plt.title( "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
    alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=15)     
#plt.legend()


non_int = hamiltonian + tunnel

values,_ = np.linalg.eig(non_int)
print "eigs:", values

print "For plotting reasons, it is wise to print these values."
print "X\t%.3f\t%.3f\t%.3f"  % (epsilon.min(), epsilon.max(), epsilon.mean())
print "Y\t%.3f\t%.3f\t%.3f" % ( transmission.min(), transmission.max(), transmission.mean())

height = values*0
if plotting_mode == 0 or plotting_mode == 2:
    height += np.average(transmission)
elif plotting_mode == 1 or plotting_mode == 3:
    height += np.average(spectral)

if decorative:
    plt.plot(values, height, 'kd', label='Eigenvalues ($H_0(V) + \\tau$)', markersize=10) 
else:
    plt.plot(values, height, 'ko', label='Eigenvalues ($H_0(V) + \\tau$)', linewidth=2) 
plt.legend()
 
plt.savefig('perspin_two_site.pdf') 

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
