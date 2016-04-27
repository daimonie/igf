import numpy as np
import scipy.interpolate as si
from scipy.optimize import minimize
from scipy.constants import physical_constants as pc
#Parallel processes!
import multiprocessing as mp
from igf import *
import sys as sys
#Command line arguments.
import argparse as argparse  
import time
#import matplotlib
#import matplotlib.pyplot as plt
from experiment import *
###

parser	= argparse.ArgumentParser(prog="rough fit",
  description = "Calculates the sum of error squared for a number of parameters..")  
  
parser.add_argument(
    '-s',
    '--sep',
    help='File number',
    action='store',
    type = int,
    default = 656
)   
parser.add_argument(
    '-c',
    '--cores',
    help='Number of cores',
    action='store',
    type = int,
    default = 4
)   
args	= parser.parse_args() 

sep = args.sep
cores = args.cores
###
epsilon_res = 1000
beta = 250
###
global_time_start = time.time()

def lsqe(x, bias_array, current_array):
    global beta, epsilon_res
    #return 1.0, 1337.666
    
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
        
        current_fit.append(xcurrent)
    ###
    current_fit = np.array(current_fit)
    scale = pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]

    current_fit *= scale
    
    scale, error,scaleerror = calculate_error( bias_array, current_fit, current_array)
    
    
    return scale, error, scaleerror
###

exp_file = "exp_data/IV130328_7_%d.dat" %sep

bias, current = read_experiment(exp_file)

#make current symmetric
current = (current - current[::-1])/2.0

points = 5 
filter = np.ones(points)/points
current = np.convolve(current, filter, mode='same')


param_list = []

for levels in np.linspace( 0.00, -0.50, 1):
    for tau in np.array([0.004]):
        for gamma in np.array([0.010]):
            for alpha in np.array([.75]):
                for capacitive in np.array([0.10]):
                    x = [tau, gamma, levels, alpha, capacitive, bias, current]
                    param_list.append(x)
                    

def error_task( argument ):
    param_x = [ argument[0], argument[1], argument[2], argument[3], argument[4] ]
    param_bias = argument[5]
    param_current = argument[6]
     
    scale, error, scaleerror =  lsqe( param_x, param_bias, param_current) 
    
    task_result = []
    task_result.extend(param_x)
    task_result.append(scale)
    task_result.append(error) 
    task_result.append(scaleerror) 
    
    #print task_result
    return task_result
                    
def super_error_task( argument_list ):
    return [ error_task(argument) for argument in argument_list]
    
parallel_pool = mp.Pool(processes=cores) #automatically uses all cores 

delta = int(len(param_list) / (cores-1))

pool_arguments = [param_list[ (delta*n):(delta*(n+1))] for n in range(cores-1)]

if len(param_list) > delta*(cores-1):
    pool_arguments.append( param_list[delta*(cores-1):len(param_list)])
 
results = parallel_pool.map(super_error_task, pool_arguments)

for worker in results:
    for row in worker:
        print "%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t" % (row[0], row[1], row[2], row[3], row[4], row[5], row[6])
###
global_time_end = time.time ()
print "\nTime spent %.6f seconds. \n " % (global_time_end - global_time_start)
