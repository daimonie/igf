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

parser  = argparse.ArgumentParser(prog="perspin transport map",
  description = "Compare spinfull/spinless transmissions.")  
  
parser.add_argument(
    '-k',
    '--ksi',
    help='Onsite interaction',
    action='store',
    type = float,
    default = 1.0
)   
parser.add_argument(
    '-u',
    '--capacitive',
    help='Onsite interaction',
    action='store',
    type = float,
    default = 0.5
)   


#change these numbers based on visual inspection
epsilon_left = -.5
epsilon_right = .5
epsilon_res = 250
epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res); 

alpha = 0.5
tau = 0.02
gamma = 0.01
bias = 0.25

args    = parser.parse_args()  
capacitive = args.capacitive
ksi = args.ksi

#capacitive = .500
#ksi = 20.00


beta = 250.0

array_levels = []
array_epsilon = []
array_perspin = []
array_perrin = []

params = np.linspace( -.5, -1e-6, 25)

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

    
    perspin_chances = perspin_calculation.distribution()
    for i in range(len(perspin_chances)):
        if perspin_chances[i] > 1e-3:
            #print i, perspin_chances[i]
            perspin_ket = perspin_calculation.ket(i)
            #print perspin_ket, perspin_ket[0]
            print "%d\t%.3f\t%d\t%d\t%d\t%d" % (i, perspin_chances[i], perspin_ket[0], perspin_ket[1], perspin_ket[2], perspin_ket[3])

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
    print "Next iteration.\n"
    ###
    print >> sys.stderr, "Finished %.3f ." % levels
    
fig, (ax1, ax2) = plt.subplots(2, figsize=(25, 15), dpi=1080)

fig.subplots_adjust(hspace=0.4)

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

max_contour = np.max([ mesh_perrin.max(), mesh_perspin.max()])

color_levels = MaxNLocator(nbins=20).tick_values(0, max_contour)
color_levels2 = MaxNLocator(nbins=20).tick_values(0, max_contour) 

cf = ax1.contourf(mesh_epsilon, mesh_param, mesh_perrin, cmap=cmap, levels=color_levels)
cf = ax2.contourf(mesh_epsilon, mesh_param, mesh_perspin, cmap=cmap, levels=color_levels2)


cb1 = fig.colorbar(cf, ax=ax1, shrink=0.9)    
cb2 = fig.colorbar(cf, ax=ax2, shrink=0.9)    

for t in cb2.ax.get_yticklabels():
     t.set_fontsize(20)
for t in cb1.ax.get_yticklabels():
     t.set_fontsize(20)
     
font = {
    'family' : 'serif',
    'weight' : 'bold',
    'size'   : 20
}

matplotlib.rc('font', **font)
#plt.rc('font', family='serif')

ax1.set_axis_bgcolor('black'); 
ax2.set_axis_bgcolor('black'); 

ax1.set_xlabel( "Energy $\\epsilon$",fontsize=30);
ax1.set_ylabel( "Zero-bias Level $\\epsilon_0$" ,fontsize=30);


ax2.set_xlabel( "Energy $\\epsilon$",fontsize=30);
ax2.set_ylabel( "Zero-level $\\epsilon_0$" ,fontsize=30);

ax1.set_title( "Spinless two-site $U=%.3f$, $\\alpha=%.3f$" % (capacitive, alpha) , fontsize=25) 
ax2.set_title( "Spinfull two-site $U=%.3f$, $\\xi=%.3f$, $V=%.3f$" % (capacitive, ksi, bias) , fontsize=25) 
 

###
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
###
plt.savefig('perspin_transport_map.pdf') 
