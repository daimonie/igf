#Numeric Python. Should be familiar.
import numpy as np
#Plotting tools
import matplotlib
import matplotlib.pyplot as plt
#Our calculation module
from igf import *
#sys. sys.stderr is a nice text stream.
import sys as sys
#Command line arguments.
import argparse as argparse  
#Time, because I want to see the total calculation time in ms.
import time
#Start the clock.
global_time_start = time.time()
print "Python version is %s.%s.%s., should be >2.7.10 for us. \n" % (sys.version_info[0],sys.version_info[1],sys.version_info[2])
#two level system, to start with
#this is the simplest couloumb-blockade system around

# Normally I use argparse for this, but due to sverely limited computational 
#   resources, i.e. a computer that should've been retired years ago, I a
#   working online (Get Data Joy) and can't use commandline arguments.   
plotting_mode = 0
chain_length = 7

capacitive_strength = .6
tunnel_strength = -0.9
epsilon_gap = .15

epsilon_left = -5.0
epsilon_right = 5.0
resolution = 100

parser	= argparse.ArgumentParser(prog="N-chain",
  description = "Calculates eigenvalues of a chain of length N.")  
  
parser.add_argument(
    '-m',
    '--mode',
    help='Plotting Mode. Zero is transmission, 1 spectral, 2 saves transmission to a svg, 3 does the same for spectral.',
    action='store',
    type = int,
    default = plotting_mode
)   
parser.add_argument(
    '-n',
    '--chain',
    help='Number of elements in the chain.',
    action='store',
    type = int,
    default = chain_length
)   
parser.add_argument(
    '-r',
    '--res',
    help='Resolution (number of horizontal data points)',
    action='store',
    type = int,
    default = resolution
)   
parser.add_argument( 
    '--cstr',
    help='Interaction Strength U',
    action='store',
    type = float,
    default = capacitive_strength
)   
parser.add_argument( 
    '--tstr',
    help='Tunnel strength (between elements)',
    action='store',
    type = float,
    default = tunnel_strength
)   
parser.add_argument( 
    '--eps',
    help='Epsilon Gap',
    action='store',
    type = float,
    default = epsilon_gap
)   
parser.add_argument( 
    '--el',
    help='Minimum energy epsilon.',
    action='store',
    type = float,
    default = epsilon_left
)   
parser.add_argument( 
    '--er',
    help='Maximum energy epsilon.',
    action='store',
    type = float,
    default = epsilon_right
)   
args	= parser.parse_args() 

plotting_mode = args.mode
chain_length = args.chain
capacitive_strength = args.cstr
tunnel_strength = args.tstr
epsilon_gap = args.eps
epsilon_left = args.el
epsilon_right = args.er
resolution = args.res

dot_levels = 2* chain_length #All at the same energy; this is a chain.

param_u = np.zeros((dot_levels, dot_levels))
param_tau = np.zeros((dot_levels, dot_levels))

for i in range(0, dot_levels):
    if i%2==0:
        l = i
        r = l + 1
        
        param_u[l][r] = capacitive_strength;
        param_u[r][l] = capacitive_strength;
        
    if i%2==0 and i/2 < chain_length-1:
        l = i
        r = i+2
        for dl in range(0,2):
            for dr in range(0,2):
                param_tau[l+dl][r+dr] = tunnel_strength;
                param_tau[r+dr][l+dl] = tunnel_strength; 
param_epsilon = np.diag( np.ones((dot_levels)))

#This makes coupling to the leads
#comparable to the coupling between levels. 
gamma_strength = 0.05 

param_gamma_left = np.zeros((dot_levels,dot_levels))
param_gamma_right = np.zeros((dot_levels,dot_levels))

param_gamma_left[0][0] = gamma_strength;
param_gamma_right[epsilon_gap][epsilon_gap] = gamma_strength;

#Inverse temperature (units same as the others)
param_beta = 0.05 * 50
 
non_int = param_epsilon + param_tau

values,_ = np.linalg.eig(non_int)

height = values*0
if plotting_mode == 0 or plotting_mode == 2:
    height += 1
elif plotting_mode == 1 or plotting_mode == 3:
    height += 1

plt.plot(values, height, 'ko') 

plt.xlim([epsilon_left, epsilon_right])
plt.ylim([0, 2])
plt.show()

global_time_end = time.time ()

print "\n Time spent %.6f seconds in mode %d. \n " % (global_time_end - global_time_start, plotting_mode)