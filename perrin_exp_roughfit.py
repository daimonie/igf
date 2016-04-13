import numpy as np
import scipy.interpolate as si
from scipy.optimize import minimize
from scipy.constants import physical_constants as pc
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
import matplotlib
import matplotlib.pyplot as plt
from experiment import *
###

epsilon_res = 1000
beta = 250
capacitive = 0.00

###
global_time_start = time.time()

def lsqe(x, bias_array, current_array):
    global beta, capacitive, epsilon_res
    
    
    xtau = x[0]
    xgamma = x[1]
    xlevels = x[2]
    xalpha = x[3]
    
    current_fit = [] 
    
    for bias in bias_array:
        epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
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

        hamiltonian[0][0] = xlevels + 0.5 * xalpha * bias
        hamiltonian[1][1] = xlevels - 0.5 * xalpha * bias

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
        
        xcurrent = np.trapz(transmission, epsilon)
        print bias, xcurrent
        
        current_fit.append(xcurrent)
    ###
    print current_fit
    scale = np.max(current) / np.max(current_fit)
    current_fit *= scale
    print current_fit
    
    
    squares = np.square( current - current_fit)
    
    sum_least_squares = squares.sum()
    
    return scale, sum_least_squares, np.array(current_fit)
###
x = np.array([0.024, 0.0102, -0.05, 0.74])

exp_file = "exp_data/IV130328_7_656.dat"

bias, current = read_experiment(exp_file)

points = 5 
filter = np.ones(points)/points
current = np.convolve(current, filter, mode='same')
        
        
scaler, error, fit_current = lsqe(x, bias, current)


print "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (x[0], x[1], x[2], x[3], scaler, error)

plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

fit_current = np.array(fit_current)
print bias.shape
print fit_current.shape
print current.shape

print fit_current

plt.plot(bias, fit_current, 'g-', )   
plt.plot(bias, current * scaler, 'r--')   

xlabel = "Bias voltage $V_B$"
ylabel = "Current $I$"
plt.rc('font', family='serif') 

plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)

plt.title( "Fitted current to experiment, scale=%.3f" % scaler)
plt.show()


global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
