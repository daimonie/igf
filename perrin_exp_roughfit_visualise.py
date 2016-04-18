import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import physical_constants as pc
from igf import *
from experiment import *
import sys as sys
import argparse as argparse  
import time
global_time_start = time.time()

plotting_mode = 0

###
sep = 648


exp_file = "exp_data/IV130328_7_%d.dat" % sep

experimental_bias, experimental_current = read_experiment(exp_file)

#make current symmetric
experimental_current = (experimental_current - experimental_current[::-1])/2.0
print "Current should be anti-symmetric; tried (I(x) - I(-x))/2 to make it look nicer.
points_filter = 5 
conv_filter = np.ones(points_filter)/points_filter
experimental_current = np.convolve(experimental_current, conv_filter, mode='same')

filename = "rough_fit_%d.txt" % sep
row = 198

file_handler = open( filename, "r" );
data = np.genfromtxt(file_handler, dtype=None,skip_footer=3);  

tau = data[row, 0]
gamma = data[row, 1]
levels = data[row, 2]
alpha = data[row, 3]
capacitive = data[row, 4]
scaler = data[row, 5]
error = data[row, 6]
### 
# left, right are now +- eV/2, see Fig 4b in Perrin(2014)
epsilon_res = 1000

bias_left = -1
bias_right = 1
bias_res = 100

interaction = np.zeros((2,2))
interaction[0][1] = capacitive
interaction[1][0] = capacitive


gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

beta = 250.00

biaswindow = np.linspace(bias_left, bias_right, bias_res)
    
current = []

for bias in biaswindow:
    hamiltonian = np.zeros((2,2))
    
    hamiltonian[0][0] = levels + 0.5 * alpha * bias
    hamiltonian[1][1] = levels - 0.5 * alpha * bias
    
    calculation = igfwl(
        hamiltonian, 
        tunnel,
        interaction, 
        gamma_left,
        gamma_right, 
        beta
    )
    epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
    
    #It is unfeasible to plot all the channels. Sum them up!
    
    transmission = calculation.full_transmission(epsilon)
 
    current.append( [ np.trapz(transmission, epsilon) ])

realscale = pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]

current = realscale * np.array(current)
minimum = 1.2 * np.min(current)
maximum = 1.2 * np.max(current)

plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

print "Experiment %.3e vs. Theory %.3e" % ( current.max(), experimental_current.max())

plt.plot(biaswindow, current, 'g-') 
plt.plot(experimental_bias, experimental_current, 'r--')   
 

title = "Current versus bias"
xlabel = "Bias $V_b$"
ylabel = "$I(V_b)$"
 
#plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

#plt.title( "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
#    alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=15)     
#plt.legend()


if plotting_mode == 2 or plotting_mode == 3:
    print "Can't save svg"
    sys.exit(0)
    plt.savefig('perrin_two_site.svg')
else:    
    plt.show()

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
