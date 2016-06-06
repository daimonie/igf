import sys as sys
#print "WARNING: CRASHES HOME COMPUTER. USE WITH CARE."
#sys.exit(0)
import numpy as np
import scipy.interpolate as si
from scipy.optimize import minimize
from scipy.constants import physical_constants as pc
from igf import *
#Command line arguments.
import argparse as argparse  
import time
from experiment import *
#griddata to format data
#3d surf plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from matplotlib.mlab import griddata 
###
global_time_start = time.time()
def kitty():
    print "Mew"
    sys.exit(0)
####
points = 5
tick_num = 50

cmap = plt.get_cmap('RdBu') 
###
separation_array = range(638, 670)
        
        
data_sep = np.array([])
data_bias = np.array([])
data_current = np.array([])

bias_set = False
bias_array = 0
for sep in separation_array:
    if sep != 645:
        file = "exp_data/IV130328_7_%d.dat" % sep
        #print "Reading [%s]" % file
        bias, current = read_experiment(file)
     

        filter = np.ones(points)/points
        current = np.convolve(current, filter, mode='same')
        
        if bias_set==False:
            bias_array = bias
            bias_set = True 
        
        current /= np.abs(current).max()
        current -= np.average(current)
        
        data_sep = np.append(data_sep, sep + 0*bias)
        data_bias = np.append(data_bias, bias)
        data_current = np.append(data_current, current)
###
fine_res = 100
fine_bias_array = np.linspace( np.min(bias_array), np.max(bias_array), fine_res) 
[mesh_bias, mesh_sep] = np.meshgrid(
    fine_bias_array,
    separation_array
)
 
mesh_current = griddata(
    data_bias,
    data_sep,
    data_current,
    mesh_bias,
    mesh_sep,
    interp='linear')
tick_min = np.min(data_current)
tick_max = np.max(data_current) 

fig, ax = plt.subplots(figsize=(25, 15), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 
     
color_levels = MaxNLocator(nbins=tick_num).tick_values(tick_min, tick_max) 

cf = ax.contourf(mesh_bias, mesh_sep, mesh_current, cmap=cmap, levels=color_levels)
cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)    

for t in cb.ax.get_yticklabels():
     t.set_fontsize(20)

plt.rc('font', family='serif')
#kitty()

ax.set_xlabel( "Bias $V$",fontsize=30); 
ax.set_ylabel( "Separation",fontsize=30); 
ax.set_title( "Filter average %d" % points,fontsize=30)

plt.savefig('perrin_experiment.pdf')

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
