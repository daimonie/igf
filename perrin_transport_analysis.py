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
level   = -0.50
beta    = 250.0
upsilon = 0.400

epsilon_left    = -0.50
epsilon_right   = 0.500
epsilon_res     = 250
epsilon_array   =  np.linspace(epsilon_left, epsilon_right, epsilon_res)

bias_left   = -1.00
bias_right  = 1.00
bias_res    = 100
bias_array = np.linspace(bias_left, bias_right, bias_res)

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
    return epsilon_array * 0
###
manager = taskManager( 4, transmission_at_params ) 
for bias in bias_array:
    param = [alpha, tau, gamma, level, beta, upsilon, epsilon_array, bias]    
    manager.add_params(param) 
###
manager.execute()
results = np.array(manager.final())

print results.shape



global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)