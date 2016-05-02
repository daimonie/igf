#Numeric Python. Should be familiar.
import numpy as np
#Plotting tools
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
from matplotlib import cm 
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
#   resources, i.e. a computer that should've been retired years ago, I am
#   working online (Get Data Joy) and can't use commandline arguments.   
plotting_mode = 2
chain_length = 3

capacitive_strength = 0.35
tunnel_strength = 0.5
epsilon_gap = -capacitive_strength/2.0# - capacitive_strength

epsilon_left = -.5
epsilon_right = 2.5
resolution = 10000

#Inverse temperature (units same as the others)
param_beta = 0.05 * 25

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
    '-b',
    '--beta',
    help='Inverse Temperature beta',
    action='store',
    type = float,
    default = param_beta
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
param_beta = args.beta
tunnel_strength = args.tstr
epsilon_gap = args.eps
epsilon_left = args.el
epsilon_right = args.er
resolution = args.res

dot_levels = chain_length #All at the same energy; this is a chain.

param_u = np.zeros((dot_levels, dot_levels))
param_tau = np.zeros((dot_levels, dot_levels))


for i in range(0, dot_levels): 
    if i + 1 < chain_length:
        l = i
        r = i+1
        for dl in range(0,1):
            for dr in range(0,1):
                param_tau[l+dl][r+dr] = tunnel_strength;
                param_tau[r+dr][l+dl] = tunnel_strength; 
                param_u[l+dl][r+dr] = capacitive_strength;
                param_u[r+dr][l+dl] = capacitive_strength; 

param_epsilon = np.diag( np.ones((dot_levels)))
#This makes coupling to the leads
#comparable to the coupling between levels. 
gamma_strength = 0.05 

param_gamma_left = np.zeros((dot_levels,dot_levels))
param_gamma_right = np.zeros((dot_levels,dot_levels))

param_gamma_left[0][0] = gamma_strength;
param_gamma_right[epsilon_gap][epsilon_gap] = gamma_strength;
 

calculation = igfwl(
    param_epsilon, 
    param_tau,
    param_u, 
    param_gamma_left,
    param_gamma_right, 
    param_beta
)
  
epsilon = np.linspace(epsilon_left, epsilon_right, resolution)

maximum = 0.02
superset = calculation.generate_superset(0)


def colour (number):
    global superset;
    return cm.afmhot( number** 2 *1.0/ 3 / len(superset))

z_indices = range(0, len(superset))
xy_vertices = []
colours = []
if plotting_mode == 0 or plotting_mode == 2:
    
    for i in superset:
        channel = calculation.transport_channel(i, epsilon)
        xy_vertices.append(list(zip(epsilon, channel))) 
        
        colours.append(colour(i))
    
        if np.max(channel) > maximum:
            maximum = np.max(channel)
    
    
    title = "Chain Transmission"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Transmission $T(\\epsilon)$"
elif plotting_mode == 1 or plotting_mode == 3:
    
    maximum = 1.2
    
    spectral = epsilon*0 
    
    for i in calculation.generate_superset(0): 
        channel = calculation.spectral_channel(i, epsilon) 
        xy_vertices.append(list(zip(epsilon, channel)))  
        colours.append(colour(i))
    
    title = "Chain Spectral"
    xlabel = "Energy $\\epsilon$"
    ylabel = "Spectral $A(\\omega)$"
    
###http://matplotlib.org/examples/mplot3d/polys3d_demo.html

fig = plt.figure(figsize=(10, 10), dpi=1080)
ax = fig.gca(projection='3d')
  
poly = PolyCollection(xy_vertices, facecolors=colours)
poly.set_alpha(0.7)
ax.add_collection3d(poly, zs=z_indices, zdir='y')
 
ax.set_xlim3d(epsilon_left, epsilon_right)
ax.set_ylim3d(0, len(superset))
ax.set_zlim3d(0, maximum)

ax.set_xlabel(xlabel, fontsize=30)
ax.set_zlabel(ylabel, fontsize=30) 

plt.title( "%d-%s: $\\beta=%.3f$, $\\epsilon_0=%.3f$, $\\Gamma=%.3f$, $\\tau=%.3f$ , $U=%.3f$" % (chain_length, title, calculation.beta,
    epsilon_gap, gamma_strength, tunnel_strength, capacitive_strength), fontsize=15)     
    
if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('spinless_chain_channels.svg')
else:    
    plt.show()
    
    
    
#time report
global_time_end = time.time ()

print "\n Time spent %.6f seconds in mode %d. \n " % (global_time_end - global_time_start, plotting_mode)