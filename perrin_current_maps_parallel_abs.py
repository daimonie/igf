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
import logging 

logging.basicConfig(level=logging.DEBUG)
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
    '-u',
    '--capacitive',
    help='Strength of interaction',
    action='store',
    type = float,
    default = 1.0
)   
args    = parser.parse_args() 

ksi = 1.0
zeta = 1.0

cores = args.cores
param_type = args.param
capacitive = args.capacitive

###
bias_left = -2.
bias_right = 2.
bias_res = 100
 

param_type = 'e'
param_left = -0.5
param_right= -0.5
param_res = 100

#param_type = 'U'
#param_left = 0
#param_right = -0.5
#param_res = 100

##param_left = -.05
##param_right = 0.05


#param_type = 'o'
#param_left = 0.00
#param_right = 5.00
#param_res = 40

tick_num = 100
tick_max = 'auto'
tick_min = 'auto'
#tick_min = -tick_max

cmap = plt.get_cmap('ocean') 
###

def ana_current(arguments): 
    #print "Non-Interacting, arguments", arguments, "\n\n"
    epsilon_res = 1000
    biaswindow  = arguments[0]
    alpha       = arguments[1]
    tau         = arguments[2]
    gamma       = arguments[3]
    capacitive  = arguments[4]
    beta        = arguments[5]
    levels      = arguments[6]
    levels_ni   = levels
    ksi         = arguments[7]
    realscale   = 2*pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]
    
    delta = lambda VV: np.sqrt( (alpha * VV )**2 + (2. * tau )**2)
    epsilon_1 = lambda VV: levels_ni - 0.5 * delta(VV)
    epsilon_2 = lambda VV: levels_ni + 0.5 * delta(VV)


    analytic_current_0 = lambda VV: realscale * gamma * ( 2. * tau)**2 / (delta(VV)**2 + gamma**2)
    analytic_current_1 = lambda VV: np.arctan( (0.5 * VV - epsilon_1(VV))/(gamma/2.0))
    analytic_current_2 = lambda VV: np.arctan( (0.5 * VV + epsilon_1(VV))/(gamma/2.0))
    analytic_current_3 = lambda VV: np.arctan( (0.5 * VV - epsilon_2(VV))/(gamma/2.0))
    analytic_current_4 = lambda VV: np.arctan( (0.5 * VV + epsilon_2(VV))/(gamma/2.0))
    analytic_current_5 = lambda VV: gamma/2.0/delta(VV) * np.log( ((0.5*VV-epsilon_1(VV))**2 + (gamma/2.0)**2)/((0.5*VV+epsilon_1(VV))**2 + (gamma/2.0)**2))
    analytic_current_6 = lambda VV: gamma/2.0/delta(VV) * np.log( ((0.5*VV-epsilon_2(VV))**2 + (gamma/2.0)**2)/((0.5*VV+epsilon_2(VV))**2 + (gamma/2.0)**2))

    analytic_current = lambda VV: analytic_current_0(VV) * ( analytic_current_1(VV)+ analytic_current_2(VV)+ analytic_current_3(VV)+ analytic_current_4(VV)+ analytic_current_5(VV)+ analytic_current_6(VV))

    current = analytic_current( biaswindow)
    #current /= np.max(current)
    return [biaswindow, current]

def calculate_current(arguments): 
    #print "Interacting, arguments", arguments, "\n\n"
    epsilon_res = 1000
    bias        = arguments[0]
    alpha       = arguments[1]
    tau         = arguments[2]
    gamma       = arguments[3]
    capacitive  = arguments[4]
    beta        = arguments[5]
    levels      = arguments[6]
    ksi         = arguments[7]
    zeta        = arguments[8]
    
    array_bias = []     
    bias_array_current = [] 
    realscale   = 2*pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]
    for bias in biaswindow:
        print "Calculating epsilon %2.3f, bias %2.3f, interaction %.3f. \n" % (levels, bias, capacitive)
        hamiltonian = np.zeros((2,2))

        hamiltonian[0][0] = levels + 0.5 * alpha * bias
        hamiltonian[1][1] = levels - 0.5 * alpha * bias

        tunnel = np.zeros((2,2))
        tunnel[0][1] = -tau
        tunnel[1][0] = -tau

        #change these numbers based on visual inspection
        epsilon_left = -0.50
        epsilon_right = 1.00
        epsilon_res = 1000

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
        epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);

        #It is unfeasible to plot all the channels. Sum them up!

        transmission = calculation.full_transmission(epsilon) 
        current = np.trapz(transmission, epsilon)
        #current = np.exp(levels*bias)
        current *= realscale 
        array_bias.append(bias) 
        bias_array_current.append(current)
        ###  
    return [array_bias, bias_array_current]
###
biaswindow = np.linspace(bias_left, bias_right, bias_res) 
param_space = np.linspace(param_left,param_right,param_res)

manager = 0
if capacitive > 0:
    manager = taskManager( cores, calculate_current ) 
else:
    manager = taskManager( cores, ana_current ) 
#print "Using analytic current..."
for param in param_space: 
    alpha = 0.75
    tau = 0.02
    gamma = 0.01
    bias = 0.25

    #capacitive = 0.40
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
    elif param_type == 'o':
        ksi = param 
    elif param_type == 'z':
        zeta = param  
    manager.add_params([biaswindow, alpha, tau, gamma, capacitive, beta, levels, ksi, zeta]) 
manager.execute()
results = manager.final()

i = 0
array_bias = []
array_param = []
array_current = []

for result in results:
    biaswindow = result[0]
    current = result[1]
    param_set = manager.get( i ) 
    param = param_set[6]
    
    if param_type == 'e':
        param = param_set[6]
    elif param_type == 'U':
        param = param_set[4]
    elif param_type == 'a':
        param = param_set[1]
    elif param_type == 'b':
        param = param_set[5]
    elif param_type == 't':
        param = param_set[2]
    elif param_type == 'g':
        param = param_set[3]
    elif param_type == 'o':
        param = param_set[7] 
        
    array_bias.extend(biaswindow)
    array_param.extend(np.array(current)*0 + param)
    array_current.extend(current)
    
    i += 1


print " %.3e < current < %.3e" % (np.min(array_current), np.max(array_current))
print " %.3e < log|current| < %.3e" % (np.min(np.log(np.abs(array_current)))/np.log( 10.00 ), np.max(np.log(np.abs(array_current)))/ np.log( 10.00 ))


array_current = np.log(np.abs( array_current ) + 1e-20) / np.log( 10.00 )

[mesh_param, mesh_bias] = np.meshgrid(
    param_space,
    biaswindow
)
 
mesh_current = griddata(
    array_param,
    array_bias,
    array_current,
    mesh_param,
    mesh_bias,
    interp='linear')

tick_max = np.max( array_current)
tick_min = np.min( array_current)

mesh_current[ mesh_current > tick_max ] = tick_max


fig, ax = plt.subplots(figsize=(25, 15), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 
     
color_levels = MaxNLocator(nbins=tick_num).tick_values(tick_min, tick_max) 

cf = ax.contourf(mesh_param, mesh_bias, mesh_current, cmap=cmap, levels=color_levels)
cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    

for t in cb.ax.get_yticklabels():
     t.set_fontsize(20)

plt.rc('font', family='serif')

ax.set_axis_bgcolor('black'); 

ax.set_ylabel( "Bias $V$" ,fontsize=30);
if param_type == 'e':
    ax.set_xlabel( "Zero-bias Level $\\epsilon_0$",fontsize=30);
elif param_type == 'U':
    ax.set_xlabel( "Capacitive Interaction Strength $U$",fontsize=30);
elif param_type == 'a':
    ax.set_xlabel( "Level-bias coupling constant $\\alpha$",fontsize=30);
elif param_type == 'b':
    ax.set_xlabel( "Inverse Temperature $\\beta$",fontsize=30);
elif param_type == 't':
    ax.set_xlabel( "Tunnel coupling $\\tau$",fontsize=30);
elif param_type == 'g':
    ax.set_xlabel( "Lead-coupling strength$\\Gamma$",fontsize=30);
elif param_type == 'o':
    ax.set_xlabel( "On-site interaction strength",fontsize=30);
ax.set_title( "$\\alpha=%.5f$, $\\tau=%.5f$, $\\Gamma=%.5f$, $\\epsilon_0=%.5f$, $\\beta=%.5f$, $U=%.5f$" % (alpha, tau, gamma, levels, beta, capacitive), fontsize=25, y=1.07) 


###
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
###

#param_type = args.param
#ksi = args.ksi
#zeta = args.zeta
plt.savefig("perrin_current_map_parallel_abs.pdf") 
