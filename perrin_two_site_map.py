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
#paper values
alpha = 0.74
tau = 0.0241
gamma = 0.0102
levels = -0.25 
capacitive = tau*0

#supplement values
alpha=0.48
tau=0.0091
gamma=0.0242
capacitive = .45*1

fermi = -4.8 #all numbers are relative to fermi (set to zero), but comparison plots are easier if you use it


tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

epsilon_right = 1#.2
epsilon_left = -1#-.7
epsilon_res = 100

bias_left = -1
bias_right = 1
bias_res = epsilon_res

cmap = plt.get_cmap('jet') 

interaction = np.zeros((2,2))
interaction[0][1] = capacitive
interaction[1][0] = capacitive

gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

beta = 12.5/5

tick_min = 0.0
tick_max = 1.0
tick_min_log = -10.00
tick_max_log = 0.00
tick_num = 100

is_log_plot = True
####################################33
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
    if is_log_plot:
        data_transmission.extend( np.log(transmission))
    else:
        data_transmission.extend( transmission )

 
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

if is_log_plot: 
    tick_min = tick_min_log
    tick_max = tick_max_log


color_levels = MaxNLocator(nbins=tick_num).tick_values(tick_min, tick_max) 

cf = ax.contourf(mesh_epsilon+fermi,mesh_bias,mesh_transmission, cmap=cmap, levels=color_levels)
cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    



for t in cb.ax.get_yticklabels():
     t.set_fontsize(20)

plt.rc('font', family='serif')

ax.set_axis_bgcolor('black'); 

ax.set_xlabel( "Energy $\\epsilon$" ,fontsize=30);
ax.set_ylabel( "Bias voltage $V$",fontsize=30);

if is_log_plot:
    ax.set_title( "Perrin two-site [log, Transmission], $\\alpha=%.5f$, $\\tau=%.5f$, $\\Gamma=%.5f$, $\\epsilon_0=%.5f$, $\\beta=%.5f$, $U=%.5f$" % (alpha, tau, gamma, levels, beta, capacitive), fontsize=25, y=1.07) 
else:
    ax.set_title( "Perrin two-site [Transmission], $\\alpha=%.5f$, $\\tau=%.5f$, $\\Gamma=%.5f$, $\\epsilon_0=%.5f$, $\\beta=%.5f$, $U=%.5f$" % (alpha, tau, gamma, levels, beta, capacitive), fontsize=25, y=1.07) 
    
plt.yticks(fontsize=30)

plt.savefig('perrin_two_site_map.png') 

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)