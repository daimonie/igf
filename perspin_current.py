import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.constants import physical_constants as pc
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
global_time_start = time.time()

plotting_mode = 1

if plotting_mode == 1:
    print "I will save the figure."
else:
    print "I will display a figure on screen."
beta        = 250.00  
gamma       = 0.010
tau         = 0.020
alpha       = 0.300
capacitive  = 0.117
alpha       = 0.25
levels      = -0.122 

parser  = argparse.ArgumentParser(prog="current map parallel",
  description = "Parallel computation of current (shape) versus some parameter.")  
  
parser.add_argument(
    '-k',
    '--ksi',
    help='Onsite interaction',
    action='store',
    type = float,
    default = 2.0
)   
parser.add_argument(
    '-z',
    '--zeta',
    help='Mutual interaction',
    action='store',
    type = float,
    default = 1.0
)   
args    = parser.parse_args()  
ksi = args.ksi
zeta = args.zeta

epsilon_left = -.3
epsilon_right = .2
epsilon_res = 1000


bias_left = -.25
bias_right = .25
bias_res = 100


biaswindow = np.linspace(bias_left, bias_right, bias_res);

spinfull_current = []
spinless_current = []
###############
for bias in biaswindow:
    print >> sys.stderr, "Calculating bias %.3f .\n" % bias
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

    
    epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
   
    spinfull_transmission = calculation.full_transmission(epsilon)  
    spinless_transmission = spinless_calculation.full_transmission(epsilon)
 
    spinfull_current.append( [ np.trapz(spinfull_transmission, epsilon) ])
    spinless_current.append( [ np.trapz(spinless_transmission, epsilon) ])
realscale = pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]
spinless_current = realscale*np.array(spinless_current)
spinfull_current = realscale*np.array(spinfull_current)

fig = None;

if plotting_mode == 1:
    fig = plt.figure(figsize=(15, 10), dpi=1080)
else:
    fig = plt.figure()
    7
ax = fig.add_subplot(111)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
fig.subplots_adjust(left=0.25)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

minimum = 1.2 * np.min([np.min(spinfull_current),np.min(spinless_current)])
maximum = 1.2 * np.max([np.max(spinfull_current),np.max(spinless_current)])
   
xlabel = "Energy $\\epsilon$"
ylabel = "Current [nA]"

nano = 1e-9
minimum /= nano
maximum /= nano


plt.plot(biaswindow, spinfull_current/nano, 'g-',label="Many-body (spinfull)")   
plt.plot(biaswindow, spinless_current/nano, 'r--',label="Many-body (spinless)")   

plt.legend(bbox_to_anchor=(0., 1.04, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "$\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$, $\\xi=%.3f$, $\\zeta=%.3f$" % (alpha, tau, gamma, levels, bias, beta, capacitive, ksi, zeta), fontsize=11)     
###
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
###
if plotting_mode == 1:
    plt.savefig('perspin_current.png')
else:    
    plt.show()
