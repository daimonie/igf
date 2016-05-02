import numpy as np
import scipy as sc
import scipy.misc as scmisc
import scipy.special as scspecial

#3d surf plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from matplotlib.mlab import griddata 

import sys as sys
#Command line arguments.
import argparse as argparse  
import time
global_time_start = time.time()

color_list = ['afmhot', 'summer','winter']
#### This is where you can make changes ;-)
electron_phonon_coupling = 1.5

n_max = 25
m_max = 25

resolutionx = int(100 * (10.00)**0.5)
resolutiony = resolutionx

x = 3
y = 1

phonon_energy = 0.015 # just a low number.
 
###end of params

def overlap(n, m, l):
    n = int(n)
    m = int(m)
    if n >= 0 and m >= 0:
        smaller =  np.min( [n, m] )
        
        fc = np.sqrt(scspecial.factorial(n) * scspecial.factorial(m))
        fc *= np.exp(-0.5 * l**2)
        
        fc *= np.sum([ (-l)**(n-k) * l**(m-k) / scspecial.factorial(k) / scspecial.factorial(m-k) / scspecial.factorial(n-k)   for k in range(smaller+1)])
        
        fc = fc
        return fc
    return 0.0
def sum_factor(n, m, l, x, y):
    global phonon_energy
    fc = overlap(n+y-x, m-y, l)
    sc = (l * phonon_energy)**x  
    sc *= scmisc.comb(x, y)
     
    sc *= fc 
    
    return sc

     
jactus = np.random.randint(0, high=len(color_list)) 
cmap = plt.get_cmap(color_list[jactus])   

number_bins = 10

[n, m] = np.meshgrid(
    np.linspace(0, n_max, resolutionx),
    np.linspace(0, m_max, resolutiony)
);
factorial = lambda xx: scmisc.factorial(xx)

z = n*m*0
#z = sum_factor(n, m, 1.1, 0, 0)
for i in range(resolutionx):
    for j in range(resolutiony):
        z[i, j] =  sum_factor(n[i, j], m[i, j], electron_phonon_coupling, x, y)

#Always make sure the final results are a meshgrid(x,y) and result z
xlabel = "Left state $n$"
ylabel = "Right state $m$"
title = "Franck-Condon, eq (6.12) of [Book Chapter], $\lambda=%.3f$, order %d, term %d, max %.3e" % (electron_phonon_coupling, x, y, z.max())

print color_list[jactus]
##don't edit here.


fig, ax = plt.subplots(figsize=(25, 15), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 


color_levels = MaxNLocator(nbins=number_bins).tick_values(z.min(), np.max([1.0, z.max()])) 

cf = ax.contourf(n,m, z, cmap=cmap, levels=color_levels)
fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.15)

plt.rc('font', family='serif')


ax.set_xlabel( "%s" % xlabel ,fontsize=30);
ax.set_ylabel( "%s" % ylabel,fontsize=30);

ax.set_title( "%s" % title, fontsize=25, y=1.07) 
    
plt.yticks(fontsize=30) 
#plt.savefig('doodling.svg')
plt.show()
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)