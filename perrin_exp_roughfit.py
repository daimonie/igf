import numpy as np
import scipy.interpolate as si
from scipy.optimize import minimize
from scipy.constants import physical_constants as pc
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
#import matplotlib
#import matplotlib.pyplot as plt
from experiment import *
###

epsilon_res = 1000
beta = 250
###
global_time_start = time.time()

def lsqe(x, bias_array, current_array):
    global beta, epsilon_res
    
    xtau = x[0]
    xgamma = x[1]
    xlevels = x[2]
    xalpha = x[3]
    capacitive = x[4]
    
    current_fit = [] 
    
    #bias_iteration = 0
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
        
        #if bias_iteration%100==0:
            #print bias_iteration, bias, xcurrent
        #bias_iteration += 1
        current_fit.append(xcurrent)
    ###
    current_fit = np.array(current_fit)
    #print current_fit
    scale = 1.00#np.max(current) / np.max(current_fit)
    current_fit *= scale
    #print current_fit
    
    
    squares = np.square( current - current_fit)
    
    sum_least_squares = squares.sum()
    
    return scale, sum_least_squares, current_fit
###

exp_file = "exp_data/IV130328_7_656.dat"

bias, current = read_experiment(exp_file)

#make current symmetric
current = (current + current[::-1])/2.0

points = 5 
filter = np.ones(points)/points
current = np.convolve(current, filter, mode='same')

fitted_current = current*0

min_error = float("Inf")
min_x = np.array([0, 0, 0, 0, 0])
min_scaler = 1.00
for levels in np.array([0.00, -0.10, -0.25, -0.45]):
    for tau in np.array([2.0, 6.0, 40.0])/1000.0:
        for gamma in np.array([10.0,100.0])/1000.0:
            for alpha in np.array([0.25, 0.50, 0.75, 1.00]):
                for capacitive in np.array([0.00, 0.10, 0.25, 0.45]):
                    x = np.array([tau, gamma, levels, alpha, capacitive])
                    scaler, error, fit_current = lsqe(x, bias, current)
                    #scaler = error = 1
                    print "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3e\t%.3e\t" % (x[0], x[1], x[2], x[3], x[4], scaler, error)
                    
                    if min_error == -1 or error < min_error:
                        print "New minimum error %.3e < %.3e" % (error, min_error)
                        fitted_current = fit_current
                        min_error = error
                        min_x = x
                        min_scaler = scaler  
###                    
#plt.figure(figsize=(10, 10), dpi=1080)
#plt.xticks(fontsize=30)
#plt.yticks(fontsize=30)
 

#plt.plot(bias, fitted_current, 'g-', )   
#plt.plot(bias, current, 'r--')   

#xlabel = "Bias voltage $V_B$"
#ylabel = "Current $I$"
#plt.rc('font', family='serif') 

#plt.xlabel(xlabel, fontsize=30)
#plt.ylabel(ylabel, fontsize=30)

#plt.title("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t" % (min_x[0], min_x[1], min_x[2], min_x[3], min_x[4], min_scaler, min_error))
#plt.savefig("roughfit.png")


global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)
