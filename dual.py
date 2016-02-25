import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
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

parser	= argparse.ArgumentParser(prog="N-chain",
  description = "Calculates tranmission or spectral function through a chain of N elements.")  
  
parser.add_argument(
    '-m',
    '--mode',
    help='Plotting Mode. Zero is transmission, 1 spectral, 2 saves transmission to a svg, 3 does the same for spectral.',
    action='store',
    type = int,
    default = plotting_mode
)   
args	= parser.parse_args() 

plotting_mode = args.mode

epsilon_gap = 0.15
gamma_strength = 0.05 
tunnel_strength = 0.00

capacitive_strength = 0.45
param_beta = 0.05 * 50


##plotting information; x lim/ylim/resolution
epsilon_left = -0.5
epsilon_right = 1.0
resolution = 500
maximum = 0.02


param_u = (np.ones(2) - np.eye(2)) * capacitive_strength
param_tau = -(np.ones(2) - np.eye(2)) * tunnel_strength
param_epsilon = np.eye(2) * epsilon_gap

#This makes coupling to the leads
#comparable to the coupling between levels. 

param_gamma_left = gamma_strength * np.eye(2)
param_gamma_right = gamma_strength * np.eye(2)

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
#Let's report probabilities

print "Temperature influences probabilities: \n"
chances = calculation.distribution ()
for i in calculation.generate_superset(0):
    print " p_%d = %.3f \n " % (i, chances[i])

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')
if plotting_mode == 0 or plotting_mode == 2:
    T0 = calculation.transport_channel(0, epsilon)
    T1 = calculation.transport_channel(1, epsilon)
    T2 = calculation.transport_channel(2, epsilon)
    T3 = calculation.transport_channel(3, epsilon) 
    
    
    plt.plot(epsilon, T0, 'r-', label="00")
    plt.plot(epsilon, T2, 'g-', label="10")
    plt.plot(epsilon, T1, 'b--', label="01")
    plt.plot(epsilon, T3, 'ko', label="11") 
    plt.plot(epsilon, T0+T1+T2+T3, 'm-', label='sum')
    
    xlabel = "energy $\\epsilon$"
    ylabel = "Transmission"
    
    if np.max(T0+T1+T2+T3) > maximum:
        maximum = np.max(T0+T1+T2+T3)*1.25
    
    plt.xlim([epsilon_left, epsilon_right])
    plt.ylim([0, maximum])
    
    title = "Single Dot, Transmission"
elif plotting_mode == 1 or plotting_mode == 3:
    A0 = calculation.spectral_channel(0, epsilon)
    A1 = calculation.spectral_channel(1, epsilon)
    A2 = calculation.spectral_channel(2, epsilon)
    A3 = calculation.spectral_channel(3, epsilon) 
    
    
    plt.plot(epsilon, A0, 'r-', label="00")
    plt.plot(epsilon, A2, 'g-', label="10")
    plt.plot(epsilon, A1, 'b--', label="01")
    plt.plot(epsilon, A3, 'ko', label="11") 
    plt.plot(epsilon, A0+A1+A2+A3, 'm-', label='sum')
    
    
    if np.max(A0+A1+A2+A3) > maximum:
        maximum = np.max(A0+A1+A2+A3)*1.25
    
    plt.xlim([epsilon_left, epsilon_right])
    plt.ylim([0, maximum])
    
    xlabel = "energy $\\epsilon$"
    ylabel = "Spectral"
     
    title = "Single Dot, Spectral" 

plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( " %s: $\\beta=%.3f$, $\\epsilon_0=%.3f$, $\\Gamma=%.3f$, $\\tau=%.3f$ , $U=%.3f$" % (title, calculation.beta,
    epsilon_gap, gamma_strength, tunnel_strength, capacitive_strength), fontsize=15)
plt.legend()
if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('dual.svg')
else:    
    plt.show()

global_time_end = time.time ()

print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)