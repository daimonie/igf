import numpy as np

 
#3d surf plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from matplotlib.mlab import griddata 

from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
global_time_start = time.time()

plotting_mode = 0

alpha = 0.74
tau = 0.0241
gamma = 0.0102
levels = -0.25 


tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

epsilon_left = - 2.0
epsilon_right = 2.0
epsilon_res = 100

bias_left = -1
bias_right = 1
bias_res = 100

interaction = np.zeros((2,2))

gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

beta = 0.05 * 50

data_bias = []
data_eps = []
data_transmission = []

epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);

for bias in np.linspace(bias_left, bias_right, bias_res):
        
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

    transmission = epsilon*0
    
    for i in calculation.generate_superset(0):
        transmission += calculation.transport_channel(i, epsilon)
    
    data_bias.extend( bias + epsilon*0)
    data_eps.extend(epsilon)
    data_transmission.extend( np.log(transmission))

 
[mesh_epsilon, mesh_bias] = np.meshgrid(
    epsilon,
    np.linspace(bias_left, bias_right, bias_res)
)
 
mesh_transmission = griddata(
    data_eps,
    data_bias,
    data_transmission,
    mesh_epsilon,
    mesh_bias,
    interp='linear')
   
fig, ax = plt.subplots(figsize=(25, 15), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 

cmap = plt.get_cmap('afmhot') 

levels = MaxNLocator(nbins=100).tick_values(mesh_transmission.min(), mesh_transmission.max()*1.1) 

cf = ax.contourf(mesh_epsilon,mesh_bias,mesh_transmission, cmap=cmap, levels=levels)
fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    

plt.rc('font', family='serif')
plt.xlim([epsilon_left, epsilon_right])

ax.set_axis_bgcolor('black'); 

ax.set_xlabel( "Energy $\\epsilon$" ,fontsize=30);
ax.set_ylabel( "Bias voltage $V$",fontsize=30);  
plt.yticks(fontsize=30)
 
plt.show()

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)