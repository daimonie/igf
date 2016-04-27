import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from igf import *
import sys as sys
#Command line arguments.
from scipy.optimize import minimize
import argparse as argparse  
import time
global_time_start = time.time()
###



separations = np.array(range(16))*2 + 638

#minrows = [41, 9, 9, 9, 43, 42, 42, 42, 42, 42, 1,1,1,1,1,1,1,1,1,1,1 ]


tau_array = []
gamma_array = []
levels_array = []
alpha_array = []
capacitive_array = []
scaler_array = []
error_array = []
scaleerror_array = []

for i in range(16):
    
    sep = separations[i]
    #row = minrows[i]
    
    filename = "rough_fit_%d.txt" % sep
    
    #print "Processing [%s]" % filename
    file_handler = open( filename, "r" );
    data = np.genfromtxt(file_handler, dtype=None,skip_footer=3);  
    
    row = np.where(data[:, 6] <= data[:, 6].min()*1.05)[0]
    row = row.min()
    
    row2 = np.where(data[:, 5] <= data[:, 5].min()*1.05)[0]
    row2 = row.min()
    
    print "%d\t%d\t%d\t%d" % (i, sep, row, row2)
    
    row = np.min([row, row2])
    
    tau = data[row, 0]
    gamma = data[row, 1]
    levels = data[row, 2]
    alpha = data[row, 3]
    capacitive = data[row, 4]
    scaler = data[row, 5]
    error = data[row, 6]
   # scaleerror = data[row, 7]

    tau_array.append(tau)
    gamma_array.append(gamma)
    levels_array.append(levels)
    alpha_array.append(alpha)
    capacitive_array.append(capacitive)
    scaler_array.append(scaler)
    error_array.append(error)
    #scaleerror_array.append(scaleerror)


###
###
plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
 
 
#title = "Scaler max %2.3e, error max %2.3e, scaleerror max %2.3e" % (np.max(scaler_array), np.max(error_array), np.max(scaleerror_array))
title = "Scaler max %2.3e, error max %2.3e" % (np.max(scaler_array), np.max(error_array))
xlabel = "Separation"
ylabel = "Value"
plt.rc('font', family='serif')


plt.plot(separations, tau_array, 'r-', label='$\\tau$')  
plt.plot(separations, gamma_array, 'g-', label='$\\Gamma$')  
plt.plot(separations, levels_array, 'b-', label='$\\epsilon_0$')  
plt.plot(separations, alpha_array, 'c-', label='$\\alpha$')  
plt.plot(separations, capacitive_array, 'm-', label='$U$')  
plt.plot(separations, scaler_array, 'k--', label='scale')  
plt.plot(separations, error_array, 'y--', label='error')  
#plt.plot(separations, scaleerror_array / np.max(scaleerror_array), 'y-', label='error')  

plt.ylim([-1.0, 1.0])
plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

plt.title( "%s" % (title), fontsize=15)     
plt.legend(loc='upper left')

global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)

plt.show() 

print "i\t t\t g\t e\t a\tU"
for i in range(16):
    print "%d\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e" % (i, tau_array[i], gamma_array[i], levels_array[i], alpha_array[i], capacitive_array[i])     
    
print "ave\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e" % (np.average(tau_array), np.average(gamma_array), np.average(levels_array), np.average(alpha_array), np.average(capacitive_array))    