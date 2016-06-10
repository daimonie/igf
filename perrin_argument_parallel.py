import numpy as np
from paralleltask import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.mlab import griddata 
from scipy.constants import physical_constants as pc
from igf import *
import sys as sys
import argparse as argparse  
###
import time
global_time_start = time.time() 
###


parser  = argparse.ArgumentParser(prog="current map parallel",
  description = "Parallel computation of current (shape) versus some parameter.")  
 
parser.add_argument(
    '-c',
    '--cores',
    help='Number of cores',
    action='store',
    type = int,
    default = 4
)   
parser.add_argument(
    '-p',
    '--param',
    help='Parameter',
    action='store',
    type = str,
    default = 'e'
)   
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
parser.add_argument(
    '-i',
    '--state',
    help='Many-body state whose eigenvalues we want to see',
    action='store',
    type = int,
    default = 0
)   
args    = parser.parse_args() 
 
cores = args.cores
param_type = args.param
ksi = args.ksi
zeta = args.zeta
state = args.state

###
bias_left = 0.00
bias_right = 1.00
bias_res = 100
 

param_type = 'e'
param_left = -1e-6
param_right = -0.5
param_res = 250


#param_type = 'o'
#param_left = 0.00
#param_right = 5.00
#param_res = 40

tick_num = 100
tick_max = 'auto'
tick_min = 'auto'
#tick_min = -tick_max

cmap = plt.get_cmap('afmhot') 
###

def calculate_argument(arguments): 
    print "Argument, arguments", arguments, "\n\n"
    epsilon_res = 250
    biaswindow  = arguments[0]
    alpha       = arguments[1]
    tau         = arguments[2]
    gamma       = arguments[3]
    capacitive  = arguments[4]
    beta        = arguments[5]
    levels      = arguments[6] 
    ksi         = arguments[7] 
    zeta        = arguments[8] 
    state_to_look_at = arguments[9]
    
    print "Look at state %d" % state_to_look_at
    angle_plus = []
    angle_minu = []
    
    for bias in biaswindow:
        
        hamiltonian = np.zeros((2,2))

        hamiltonian[0][0] = levels + 0.5 * alpha * bias
        hamiltonian[1][1] = levels - 0.5 * alpha * bias

        tunnel = np.zeros((2,2))
        tunnel[0][1] = -tau
        tunnel[1][0] = -tau

        #change these numbers based on visual inspection
        epsilon_left = -0.50
        epsilon_right = 1.00
        epsilon_res = 10000

        interaction = np.zeros((2,2))
        interaction[0][1] = capacitive
        interaction[1][0] = capacitive


        gamma_left = np.zeros((2,2))
        gamma_left[0][0] = gamma

        gamma_right = np.zeros((2,2))
        gamma_right[1][1] = gamma

        beta = 250.0*1

        calculation = igfwl(
            hamiltonian, 
            tunnel,
            interaction, 
            gamma_left,
            gamma_right, 
            beta
        )
        
        state = calculation.ket(state_to_look_at)
        capacitive_contribution = np.dot(state.T, np.dot( interaction, state))
        
        angle_plus.append(levels + capacitive_contribution + 0.5*alpha * bias + 1j*gamma)
        angle_minu.append(levels + capacitive_contribution - 0.5*alpha * bias + 1j*gamma)    
        
    angle_plus = np.angle(angle_plus)
    angle_minu = np.angle(angle_minu)
    return [biaswindow, angle_plus, angle_minu]

biaswindow = np.linspace(bias_left, bias_right, bias_res) 
param_space = np.linspace(param_left,param_right,param_res)

manager = taskManager( cores, calculate_argument ) 
#manager = taskManager( cores, ana_current ) 
#print "Using analytic current..."
for param in param_space: 
    alpha = 0.75
    tau = 0.02
    gamma = 0.01
    levels = -0.5
    capacitive = 0.40
    beta = 250.00
     
    levels =  -1e-6
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
    elif param_type == 'o':
        ksi = param 
    elif param_type == 'z':
        zeta = param  
    manager.add_params([biaswindow, alpha, tau, gamma, capacitive, beta, levels, ksi, zeta, state]) 
manager.execute()

print "Calculation is performed, working on countours"

results = np.array(manager.final())
angle_plus = results[:, 0]
angle_minu = results[:, 1]

[mesh_bias, mesh_param] = np.meshgrid(
    biaswindow,
    param_space
)
 
mesh_plus = angle_plus * np.pi
mesh_minu = angle_minu * np.pi

fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(20, 23), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 
     
colors_plus = MaxNLocator(nbins=tick_num).tick_values(-1.0*np.pi, 1.0*np.pi) 
colors_minu = MaxNLocator(nbins=tick_num).tick_values(-1.0*np.pi, 1.0*np.pi) 

cf_plus = ax1.contourf(mesh_bias, mesh_param, mesh_plus, cmap=cmap, levels=colors_plus)
cf_minu = ax2.contourf(mesh_bias, mesh_param, mesh_minu, cmap=cmap, levels=colors_minu)


cb_plus = fig.colorbar(cf_plus, ax=ax1, shrink=0.9, pad=0.15)    
cb_minu = fig.colorbar(cf_minu, ax=ax2, shrink=0.9, pad=0.15)    

plt.rc('font', family='serif')

for t in cb_plus.ax.get_yticklabels():
     t.set_fontsize(20)
for t in cb_minu.ax.get_yticklabels():
     t.set_fontsize(20)

ax2.set_xlabel( "Bias $V$" ,fontsize=30); 
ax1.set_ylabel( "Zero-bias Level $\\epsilon_0$",fontsize=30); 
ax2.set_ylabel( "Zero-bias Level $\\epsilon_0$",fontsize=30);

ax1.set_title( "$\\alpha=%.5f$, $\\tau=%.5f$, $\\Gamma=%.5f$, $\\epsilon_0=%.5f$, $\\beta=%.5f$, $U=%.5f$" % (alpha, tau, gamma, levels, beta, capacitive), fontsize=25, y=1.07) 

###
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
###

#param_type = args.param
#ksi = args.ksi
#zeta = args.zeta
plt.savefig("perrin_argument_parallel.pdf") 



