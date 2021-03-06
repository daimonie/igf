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

clip_size = 25.00
epsilon_gap = 0.35
gamma_strength = 0.05 
tunnel_strength = 0.00
 
param_beta = 0.05 * 50

matrix_dimensions = 4


##plotting information; x lim/ylim/resolution
epsilon_left = -0.5
epsilon_right = 3.0
resolution = 50
maximum = 0.02


param_tau = -(np.ones(matrix_dimensions) - np.eye(matrix_dimensions)) * tunnel_strength
param_epsilon = np.eye(matrix_dimensions) * epsilon_gap

print "Assume all states have spacing epsilon_gap, starting from zero. \n"

i = 0
for j in range(0, matrix_dimensions):
    if j % 2 == 0:
        i += 1
    
    param_epsilon[j][j] += epsilon_gap*(i-1)


#This makes coupling to the leads
#comparable to the coupling between levels. 

param_gamma_left = gamma_strength * np.eye(matrix_dimensions)
param_gamma_right = gamma_strength * np.eye(matrix_dimensions)


data_epsilon = []
data_capacitive = []
data_z = [] 

for capacitive_strength in np.linspace(0,1,resolution):
    param_u = (np.ones(matrix_dimensions) - np.eye(matrix_dimensions)) * capacitive_strength

    
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
       
        data_epsilon.extend(epsilon)
        data_capacitive.extend(capacitive_strength + epsilon*0)
        data_z.extend(transmission) 
    elif plotting_mode == 1 or plotting_mode == 3:
        
        
        spectral = calculation.spectral_channel(0, epsilon)
        
        for i in calculation.generate_superset(1):
            spectral += calculation.spectral_channel(i, epsilon)
         
        
        data_epsilon.extend(epsilon)
        data_capacitive.extend(capacitive_strength + epsilon*0)
        data_z.extend(spectral)
    
 #continue plotting
[mesh_epsilon, mesh_capacitive] = np.meshgrid(
    np.linspace(epsilon_left,epsilon_right, resolution),
    np.linspace(0,1,resolution)
)
 
mesh_transmission = griddata(
    data_epsilon,
    data_capacitive,
    data_z,
    mesh_epsilon,
    mesh_capacitive,
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

cf = ax.contourf(mesh_epsilon,mesh_capacitive,mesh_transmission, cmap=cmap, levels=levels)
fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    

plt.rc('font', family='serif')
plt.xlim([epsilon_left, epsilon_right])

ax.set_axis_bgcolor('black'); 

ax.set_xlabel( "Energy $\\epsilon$" ,fontsize=30);
ax.set_ylabel( "Capacitive Interaction Strength $U$",fontsize=30); 

title = ""
if plotting_mode == 0 or plotting_mode == 2:
    title = "WBL Transmission, %d-dot" % int(matrix_dimensions/2)
elif plotting_mode == 1:
    title = "WBL Spectral" 

plt.title( " %s: $\\beta=%.3f$, $\\epsilon_0=%.3f$, $\\Gamma=%.3f$, $\\tau=%.3f$" % (title, calculation.beta,
    epsilon_gap, gamma_strength, tunnel_strength), fontsize=15)

if plotting_mode == 2 or plotting_mode == 3:
    filename = 'quad_3d.svg'
    print "Saving %s ." % filename
    plt.savefig( filename )
else:    
    plt.show()


global_time_end = time.time ()

print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)