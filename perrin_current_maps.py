import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.mlab import griddata 
from igf import *
import sys as sys
import argparse as argparse  
###
import time
global_time_start = time.time() 
###

bias_left = -1
bias_right = 1
bias_res = 100

epsilon_res = 250

param_type = 'e'
param_left = -1e-6
param_right = -0.5
param_res = 40


array_param = []
array_bias = []    
array_current = [] 

tick_num = 100
tick_max = 'auto'
tick_min = 'auto'
#tick_min = -tick_max

cmap = plt.get_cmap('afmhot') 
###
biaswindow = np.linspace(bias_left, bias_right, bias_res) 
param_space = np.linspace(param_left,param_right,param_res)

for param in param_space:
    print "Param[%s]:\t%.3f" % (param_type, param)
    
    bias_array_current = [] 
    for bias in biaswindow:

        alpha = 0.74
        tau = 0.0241
        gamma = 0.0102
        levels = -0.05
        capacitive = 0.15*0
        beta = 250.00
        
        if param_type == 'e':
            levels = param
        elif param_type == 'U':
            capacitive = param
        elif param_type == 'a':
            alpha = param
        elif param_type == 'b':
            beta = param
        elif param_type == 't':
            tau = param
        elif param_type == 'g':
            gamma = param
        
        tunnel = np.zeros((2,2))
        tunnel[0][1] = -tau
        tunnel[1][0] = -tau 

        interaction = np.zeros((2,2))
        interaction[0][1] = capacitive
        interaction[1][0] = capacitive


        gamma_left = np.zeros((2,2))
        gamma_left[0][0] = gamma

        gamma_right = np.zeros((2,2))
        gamma_right[1][1] = gamma


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
        
        current = np.trapz(transmission, epsilon)
        
        array_bias.append(bias)
        array_param.append(param)
        bias_array_current.append(current)
        ###
        
    bias_array_current = np.array(bias_array_current) / np.max(bias_array_current)
    
    array_current.extend(bias_array_current)
         

[mesh_bias, mesh_param] = np.meshgrid(
    biaswindow,
    param_space
)
 
mesh_current = griddata(
    array_bias,
    array_param,
    array_current,
    mesh_bias,
    mesh_param,
    interp='linear')

if(tick_max == 'auto'): 
    tick_max = np.max(mesh_current)
    tick_min = np.min(mesh_current)
fig, ax = plt.subplots(figsize=(25, 15), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 
     
color_levels = MaxNLocator(nbins=tick_num).tick_values(tick_min, tick_max) 

cf = ax.contourf(mesh_bias, mesh_param, mesh_current, cmap=cmap, levels=color_levels)
cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    

for t in cb.ax.get_yticklabels():
     t.set_fontsize(20)

plt.rc('font', family='serif')

ax.set_axis_bgcolor('black'); 

ax.set_xlabel( "Bias $V$" ,fontsize=30);
if param_type == 'e':
    ax.set_ylabel( "Zero-bias Level $\\epsilon_0$",fontsize=30);
elif param_type == 'U':
    ax.set_ylabel( "Capacitive Interaction Strength $U$",fontsize=30);
elif param_type == 'a':
    ax.set_ylabel( "Level-bias coupling constant $\\alpha$",fontsize=30);
elif param_type == 'b':
    ax.set_ylabel( "Inverse Temperature $\\beta$",fontsize=30);
elif param_type == 't':
    ax.set_ylabel( "Tunnel coupling $\\tau$",fontsize=30);
elif param_type == 'g':
    ax.set_ylabel( "Lead-coupling strength$\\Gamma$",fontsize=30);
ax.set_title( "$\\alpha=%.5f$, $\\tau=%.5f$, $\\Gamma=%.5f$, $\\epsilon_0=%.5f$, $\\beta=%.5f$, $U=%.5f$" % (alpha, tau, gamma, levels, beta, capacitive), fontsize=25, y=1.07) 


###
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
###

plt.savefig('perrin_current_map.png') 