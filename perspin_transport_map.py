import numpy as np
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

plotting_mode = 0

#change these numbers based on visual inspection
epsilon_left = -.6
epsilon_right = .4
epsilon_res = 10000
epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res); 

alpha = 0.25
tau = 0.02
gamma = 0.01
bias = 0.0
capacitive = .117
ksi = 3.50
beta = 250.0

array_levels = []
array_epsilon = []
array_perspin = []
array_perrin = []

params = np.linspace( -.75, .5, int(1+40./2.*1.25))

for levels in params:
    print >> sys.stderr, "Starting %.3f ." % levels
    perspin_hamiltonian = np.zeros((4,4))

    perspin_hamiltonian[0][0] = levels + 0.5 * alpha * bias
    perspin_hamiltonian[1][1] = levels + 0.5 * alpha * bias
    perspin_hamiltonian[2][2] = levels - 0.5 * alpha * bias
    perspin_hamiltonian[3][3] = levels - 0.5 * alpha * bias

    perspin_tunnel = np.zeros((4,4))
    
    #left-right 
    perspin_tunnel[0][2] = tau  
    perspin_tunnel[1][3] = tau  
    
    #right-left
    perspin_tunnel[2][0] = tau  
    perspin_tunnel[3][1] = tau  

    perspin_interaction = np.zeros((4,4))

    perspin_interaction[0][1] = ksi * capacitive
    perspin_interaction[1][0] = ksi * capacitive

    perspin_interaction[2][3] = ksi * capacitive
    perspin_interaction[3][2] = ksi * capacitive

    perspin_interaction[0][3] = capacitive
    perspin_interaction[3][0] = capacitive

    perspin_interaction[0][2] = capacitive
    perspin_interaction[2][0] = capacitive

    perspin_interaction[1][3] = capacitive
    perspin_interaction[3][1] = capacitive

    perspin_interaction[1][2] = capacitive
    perspin_interaction[2][1] = capacitive
 
    perspin_gamma_left = np.zeros((4,4))
    perspin_gamma_left[0][0] = gamma
    perspin_gamma_left[1][1] = gamma

    perspin_gamma_right = np.zeros((4,4))
    perspin_gamma_right[2][2] = gamma
    perspin_gamma_right[3][3] = gamma
 
    perspin_calculation = igfwl(
        perspin_hamiltonian, 
        perspin_tunnel,
        perspin_interaction, 
        perspin_gamma_left,
        perspin_gamma_right, 
        beta
    ) 

    perspin_transmission = np.abs(perspin_calculation.full_transmission(epsilon))

    perrin_hamiltonian = np.zeros((2,2))

    perrin_hamiltonian[0][0] = levels + 0.5 * alpha * bias
    perrin_hamiltonian[1][1] = levels - 0.5 * alpha * bias

    perrin_tunnel = np.zeros((2,2))
    perrin_tunnel[0][1] = -tau
    perrin_tunnel[1][0] = -tau

    perrin_interaction = np.zeros((2,2))
    perrin_interaction[0][1] = capacitive
    perrin_interaction[1][0] = capacitive


    perrin_gamma_left = np.zeros((2,2))
    perrin_gamma_left[0][0] = gamma

    perrin_gamma_right = np.zeros((2,2))
    perrin_gamma_right[1][1] = gamma

    perrin_calculation = igfwl(
        perrin_hamiltonian, 
        perrin_tunnel,
        perrin_interaction, 
        perrin_gamma_left,
        perrin_gamma_right, 
        beta
    )
    perrin_transmission = np.abs(perrin_calculation.full_transmission(epsilon)) 
    
    array_epsilon.extend(epsilon + 0*perrin_transmission)
    array_levels.extend(levels + 0*perrin_transmission)
    
    array_perrin.extend(perrin_transmission)
    array_perspin.extend(perspin_transmission)
    
    ###chance inspection spinless
    p = perrin_calculation.distribution()
    max_p = 0.00
    max_x = perrin_calculation.ket(0)
    max_i = 0
    
    for i in perrin_calculation.generate_superset(0):
        x = perrin_calculation.ket(i)
        
        if p[i] > max_p:
            max_p = p[i]
            max_x = x
            max_i = i
        #
    #
    print "%d\t%.3e\t%d\t%d" % (max_i, max_p, max_x[0], max_x[1])
    
    ###chance inspection spin
    p = perspin_calculation.distribution()
    max_p = 0.00
    max_x = perspin_calculation.ket(0)
    max_i = 0
    
    for i in perspin_calculation.generate_superset(0):
        x = perspin_calculation.ket(i)
        
        if p[i] > max_p:
            max_p = p[i]
            max_x = x
            max_i = i
        #
    #
    print "%d\t%.3e\t%d\t%d\t%d\t%d" % (max_i, max_p, max_x[0], max_x[1], max_x[2], max_x[3])
    ###
    print >> sys.stderr, "Finished %.3f ." % levels
    
fig, (ax1, ax2) = plt.subplots(2, figsize=(25, 15), dpi=1080)

fig.subplots_adjust(hspace=0.2)

plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 


[mesh_epsilon, mesh_param] = np.meshgrid(
    epsilon, 
    params
)     

mesh_perrin = griddata(
    array_epsilon,
    array_levels,
    array_perrin,
    mesh_epsilon,
    mesh_param,
    interp='linear')
mesh_perspin = griddata(
    array_epsilon,
    array_levels,
    array_perspin,
    mesh_epsilon,
    mesh_param,
    interp='linear')
###
# output, not plot, because errors.


cmap = plt.get_cmap('afmhot') 

color_levels = MaxNLocator(nbins=20).tick_values(0, np.max([ mesh_perrin.max(), mesh_perspin.max()])) 
cf = ax1.contourf(mesh_epsilon, mesh_param, mesh_perrin, cmap=cmap, levels=color_levels)
cf = ax2.contourf(mesh_epsilon, mesh_param, mesh_perspin, cmap=cmap, levels=color_levels)


cb1 = fig.colorbar(cf, ax=ax1, shrink=0.9)    
cb2 = fig.colorbar(cf, ax=ax2, shrink=0.9)    

for t in cb2.ax.get_yticklabels():
     t.set_fontsize(20)
for t in cb1.ax.get_yticklabels():
     t.set_fontsize(20)
     
     
plt.rc('font', family='serif')

ax1.set_axis_bgcolor('black'); 
ax2.set_axis_bgcolor('black'); 

ax1.set_xlabel( "Energy $\\epsilon$",fontsize=30);
ax1.set_ylabel( "Zero-bias Level $\\epsilon_0$" ,fontsize=30);


ax2.set_xlabel( "Energy $\\epsilon$",fontsize=30);
ax2.set_ylabel( "Zero-level $\\epsilon_0$" ,fontsize=30);

ax1.set_title( "Spinless two-site $U=%.3f$" % (capacitive) , fontsize=25) 
ax2.set_title( "Spinfull two-site $U=%.3f$, $\\xi=%.3f$" % (capacitive, ksi) , fontsize=25) 
 

###
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
###
plt.savefig('perspin_transport_map.pdf') 
