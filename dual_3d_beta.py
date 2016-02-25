import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
#3d surf plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from matplotlib.mlab import griddata 
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

clip_size = 0.15

epsilon_gap = 0.15
gamma_strength = 0.05 
tunnel_strength = 0.00
  
capacitive_strength = 0.45

##plotting information; x lim/ylim/resolution
epsilon_left = -0.5
epsilon_right = 1.0
resolution = 50
maximum = 0.02

 
param_u = (np.ones(2) - np.eye(2)) * capacitive_strength
param_tau = -(np.ones(2) - np.eye(2)) * tunnel_strength
param_epsilon = np.eye(2) * epsilon_gap

max_beta = 100.0
#This makes coupling to the leads
#comparable to the coupling between levels. 

param_gamma_left = gamma_strength * np.eye(2)
param_gamma_right = gamma_strength * np.eye(2)



data_epsilon = []
data_beta = []
data_z = [] 

for param_beta in np.linspace(0,max_beta,resolution):

    
    calculation = igfwl(
        param_epsilon, 
        param_tau,
        param_u, 
        param_gamma_left,
        param_gamma_right, 
        param_beta
    )
    
    epsilon = np.linspace(epsilon_left, epsilon_right,resolution)
    if plotting_mode == 0 or plotting_mode == 2:
        
        #It is unfeasible to plot all the channels. Sum them up!
        
        transmission = calculation.transport_channel(0, epsilon)
        
        for i in calculation.generate_superset(1):
            transmission += calculation.transport_channel(i, epsilon)
       
        data_z.extend(transmission) 
    elif plotting_mode == 1 or plotting_mode == 3:
        
        
        spectral = calculation.spectral_channel(0, epsilon)
        
        for i in calculation.generate_superset(1):
            spectral += calculation.spectral_channel(i, epsilon)
         
         
        data_z.extend(spectral)
    data_epsilon.extend(epsilon)
    data_beta.extend(param_beta + epsilon*0)
    
 #continue plotting
[mesh_epsilon, mesh_beta] = np.meshgrid(
    np.linspace(epsilon_left,epsilon_right, resolution),
    np.linspace(0,max_beta,resolution)
)
 
mesh_transmission = griddata(
    data_epsilon,
    data_beta,
    data_z,
    mesh_epsilon,
    mesh_beta,
    interp='linear')
    
    
fig, ax = plt.subplots(figsize=(25, 15), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

if plotting_mode == 0 or plotting_mode == 2:
    mesh_transmission = np.clip(mesh_transmission, 0, clip_size)
elif plotting_mode == 1 or plotting_mode == 3:
    mesh_transmission = np.clip(mesh_transmission, 0, clip_size)

cmap = plt.get_cmap('afmhot') 

levels = MaxNLocator(nbins=100).tick_values(mesh_transmission.min(), mesh_transmission.max()*1.1) 

cf = ax.contourf(mesh_epsilon,mesh_beta,mesh_transmission, cmap=cmap, levels=levels)
fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    

plt.rc('font', family='serif')
plt.xlim([epsilon_left, epsilon_right])

ax.set_axis_bgcolor('black'); 

ax.set_xlabel( "Energy $\\epsilon$" ,fontsize=30);
ax.set_ylabel( "Inverse Temperature $\\beta$",fontsize=30); 

title = ""
if plotting_mode == 0 or plotting_mode == 2:
    title = "WBL Transmission, Dual, "
elif plotting_mode == 1 or plotting_mode == 3:
    title = "WBL Spectral, Dual, " 

plt.title( " %s: $\\beta=%.3f$, $\\epsilon_0=%.3f$, $\\Gamma=%.3f$, $\\tau=%.3f$, $U=%.3f$" % (title, calculation.beta,
    epsilon_gap, gamma_strength, tunnel_strength, capacitive_strength), fontsize=15)

if plotting_mode == 2 or plotting_mode == 3:
    plt.savefig('dual_3d_beta.svg')
else:    
    plt.show()


global_time_end = time.time ()

print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)