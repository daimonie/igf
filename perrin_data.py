import numpy as np
import scipy.interpolate as si
from scipy.optimize import minimize
from scipy.constants import physical_constants as pc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata 
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
global_time_start = time.time()

def read_experiment_data( filename ):
    
    file_handler = open( filename, "r" );
    
    data = np.genfromtxt(file_handler, skip_header=9, dtype=None, usecols=range(0,3)); #excluding the symtype col
    exp_bias = data[:,0]
    exp_current = data[:,1] 
    
    truebias = np.linspace(np.min(exp_bias), np.max(exp_bias), int(exp_bias.shape[0]/2))
    
    truecurrent = si.griddata(exp_bias, exp_current, truebias, method='nearest')
    
    return truebias, truecurrent
def read_experiment( filename ):
    exp_bias, exp_current = read_experiment_data(first_file)
    
    dI = np.max(exp_current) - np.min(exp_current)
    dV = np.max(exp_bias) - np.min(exp_bias)
    
    exp_background = lambda V: dI/dV * V
    
    exp_current -= exp_background(exp_bias)
    
    return exp_bias, exp_current/1e-9

first_file = "exp_data/IV130328_7_656.dat"

bias, current = read_experiment(first_file)

#these are fixed
beta = 250.00 
capacitive = 0.15
alpha = 0.74
tau = 0.0241
gamma = 0.0102
levels = -0.05
epsilon_res = 1000
#fitting time
def lsqe(x, bias_array, current_array):
    global beta, capacitive, alpha, epsilon_res
    
    xtau = x[0]
    xgamma = x[1]
    xlevels = x[2]
    xscale = x[3] * pc["elementary charge over h"][0]
    current_fit = []
    print x
    for bias in bias_array:
        epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
	print "\t V=%.3f" % bias        
        tunnel = np.zeros((2,2))
        tunnel[0][1] = -xtau
        tunnel[1][0] = -xtau 

        interaction = np.zeros((2,2))
        interaction[0][1] = capacitive
        interaction[1][0] = capacitive


        gamma_left = np.zeros((2,2))
        gamma_left[0][0] = xgamma

        gamma_right = np.zeros((2,2))
        gamma_right[1][1] = xgamma


        hamiltonian = np.zeros((2,2))

        hamiltonian[0][0] = xlevels + 0.5 * alpha * bias
        hamiltonian[1][1] = xlevels - 0.5 * alpha * bias

        calculation = igfwl(
        hamiltonian, 
        tunnel,
        interaction, 
        gamma_left,
        gamma_right, 
        beta
        )
       
        #It is unfeasible to plot all the channels. Sum them up!

        transmission = calculation.full_transmission(epsilon)
        
        current = xscale*np.trapz(transmission, epsilon)
        current_fit.append([current])
    ###
    
    squares = np.square( current - current_fit)
    
    sum_least_squares = squares.sum()
    
    return sum_least_squares
###
optimizer = minimize( lsqe, [tau, gamma, levels, 1.0], args=(bias, current),
    jac=False, method='L-BFGS-B', options={'disp': False, 'eps': np.array([tau/20.00, gamma/20.00, levels, 1e-5])}, tol = 1e-6)
    
x = optimizer.x
print "Finally: ", x  

plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

xlabel = "Bias $V$"
ylabel = "Current $\\left[nA\\right]$"
title  = "%s, removed linear background" %  first_file
plt.rc('font', family='serif')
 
 
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "%s" % (title), fontsize=15)     
 
plt.plot(bias, current, 'g-', label='experiment') 


plt.savefig("perrin_data.png")
    
###############
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
