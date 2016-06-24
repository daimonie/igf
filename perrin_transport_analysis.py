import numpy as np
from paralleltask import *
import matplotlib
#matplotlib.use('Agg')
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
#Command line arguments.
import argparse as argparse  
import time
global_time_start = time.time()

alpha   = 0.750
tau     = 0.020
gamma   = 0.010
level   = -0.2
beta    = 250.0
upsilon = 0.400

epsilon_left    = -0.50
epsilon_right   = 0.500
epsilon_res     = 250
epsilon_array   =  np.linspace(epsilon_left, epsilon_right, epsilon_res)

bias_left   = -0.00
bias_right  = 1.0
bias_res    = 100
bias_array = np.linspace(bias_left, bias_right, bias_res)

cutoff = 0.01
cmap = plt.get_cmap('afmhot') 
### parallel function
def transmission_at_params( params ):
    alpha   = params[0]
    tau     = params[1]
    gamma   = params[2]
    level   = params[3]
    beta    = params[4]
    upsilon = params[5]
    epsilon_array = params[6] #might seem odd. But we are calculating this in parallel, so I want to avoid shared memory.
    
    bias    = params[7]
    print "Calculating bias V=%2.3f" % bias
        
        ###make calculation object
        
    hamiltonian = np.zeros((2,2))

    hamiltonian[0][0] = level + 0.5 * alpha * bias
    hamiltonian[1][1] = level - 0.5 * alpha * bias

    tunnel = np.zeros((2,2))
    tunnel[0][1] = -tau
    tunnel[1][0] = -tau
 
    interaction = np.zeros((2,2))
    interaction[0][1] = upsilon
    interaction[1][0] = upsilon


    gamma_left = np.zeros((2,2))
    gamma_left[0][0] = gamma

    gamma_right = np.zeros((2,2))
    gamma_right[1][1] = gamma

    calculation = igfwl(
        hamiltonian, 
        tunnel,
        interaction, 
        gamma_left,
        gamma_right, 
        beta
    )
    ###calculate 
    
    transmission = calculation.full_transmission(epsilon_array)
    
    return transmission
### parallel task distribution
manager = taskManager( 4, transmission_at_params ) 
for bias in bias_array:
    param = [alpha, tau, gamma, level, beta, upsilon, epsilon_array, bias]    
    manager.add_params(param) 
###
manager.execute()
results = np.array(manager.final())

[mesh_epsilon, mesh_bias] = np.meshgrid(
    epsilon_array,
    bias_array
)

grid_epsilon        = []
grid_bias           = []
grid_transmission   = []

for row in range(int(results.shape[0])):
    bias = bias_array[row]
    transmission = results[row]
 
    tmax = cutoff * np.max(transmission)
    
    transmission[transmission>tmax]=tmax
    
    grid_epsilon.extend( epsilon_array )
    grid_bias.extend( bias + epsilon_array*0)
    grid_transmission.extend( transmission )

mesh_transmission = griddata(
    grid_epsilon,
    grid_bias,
    grid_transmission,
    mesh_epsilon,
    mesh_bias,
    interp='linear')
###plot time
fig, ax = plt.subplots(figsize=(25, 15), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 
     
color_levels = MaxNLocator(nbins=100).tick_values(0.00, np.max( grid_transmission )) 

cf = ax.contourf(mesh_epsilon, mesh_bias, mesh_transmission, cmap=cmap, levels=color_levels)
cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    

plt.plot( bias_array/2.0, bias_array, 'b-')   
plt.plot( -bias_array/2.0, bias_array, 'g-')   

for t in cb.ax.get_yticklabels():
     t.set_fontsize(20)

plt.rc('font', family='serif')
plt.xlim( [epsilon_left, epsilon_right])
plt.ylim( [bias_left, bias_right])

plt.xticks(np.array(range(11))*0.1-0.5)
plt.yticks(np.array(range(6))*0.2)


ax.set_axis_bgcolor('black'); 

ax.set_xlabel( "Energy $\\epsilon$" ,fontsize=30);
ax.set_ylabel( "Bias voltage $V$",fontsize=30);
ax.set_title( "Transport Analysis, $\\epsilon_0=%.3f,U=%.3f$, cutoff at %d percent, integration window green to blue." % (level, upsilon, int(cutoff*100)), fontsize=25, y=1.07) 

###elapsed time, save

plt.savefig('perrin_transport_analysis.pdf') 
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
