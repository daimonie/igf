import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
from scipy.constants import physical_constants as pc
from igf import *
from experiment import *
import sys as sys
import argparse as argparse  
import time

#Parallel processes!
import multiprocessing as mp

global_time_start = time.time()

plotting_mode = 0

###
sep = 650
row = -60
calculate_current = True
#calculate_current = False

editing_parameters = False
editing_parameters = True

###
exp_file = "exp_data/IV130328_7_%d.dat" % sep
print "Data from [%s], anti-symmetrised for a better fit." % exp_file
experimental_bias, experimental_current = read_experiment(exp_file)
### 
#make current symmetric
original_current = experimental_current
experimental_current = (experimental_current - experimental_current[::-1])/2.0

points_filter = 5 
conv_filter = np.ones(points_filter)/points_filter
experimental_current = np.convolve(experimental_current, conv_filter, mode='same')

filename = "fit/rough_fit_%d.txt" % sep

file_handler = open( filename, "r" );
data = np.genfromtxt(file_handler, dtype=None,skip_footer=3);  
if row < 0:
    error_array = data[:, 7]
    row = np.where( error_array <= error_array.min()*1.05)[0]
    row = row.min()

    print "Lowest error row is %d." % row
else:
    print "Calculate error row %d." %  row
tau = data[row, 0]
gamma = data[row, 1]
levels = data[row, 2]
alpha = data[row, 3]
capacitive = data[row, 4]
scaler = data[row, 5]
error = data[row, 6]

print "Fitted parameters:"
print "\t tau = %2.3f" % tau
print "\t gamma = %2.3f" % gamma
print "\t levels = %2.3f" % levels
print "\t alpha = %2.3f" % alpha
print "\t capacitive = %2.3f" % capacitive
print "\t scaler = %2.3f" % scaler
print "\t error = %2.3f" % error
### 
# left, right are now +- eV/2, see Fig 4b in Perrin(2014)
epsilon_res = 1000


levels_ni = 0.00
if editing_parameters:
    print "Changing parameters for my convenience."
     
    gamma = 0.010
    tau = 0.020
    alpha = 0.300
    capacitive = 0.117
    levels = -0.122
     
    alpha = 0.25
     
    levels_ni = 0.00
    
bias_left = -1
bias_right = 1
bias_res = 100

interaction = np.zeros((2,2))
interaction[0][1] = capacitive
interaction[1][0] = capacitive


gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

beta = 250.00

biaswindow = experimental_bias
realscale = pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]
print "Using %d cores" % mp.cpu_count()
current = []
if calculate_current:
    #for bias in biaswindow:
    def current_task(bias):
        global realscale, levels, alpha, tunnel, interaction, gamma_left, gamma_right, beta, epsilon_res
        hamiltonian = np.zeros((2,2))
        
        hamiltonian[0][0] = levels + 0.5 * alpha * bias
        hamiltonian[1][1] = levels - 0.5 * alpha * bias
        
        calculation = igfwl(
            hamiltonian, 
            tunnel,
            interaction, 
            gamma_left,
            gamma_right, 
            beta
        )
        epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
        
        #It is unfeasible to plot all the channels. Sum them up!
        
        transmission = calculation.full_transmission(epsilon)
        this_current = np.trapz(transmission, epsilon) 
        this_current = realscale * np.array(this_current)
        return this_current
    #    current.append( [ np.trapz(transmission, epsilon) ])
    def super_current_task(arg):
        small_window = arg[0]
        n = arg[1]
        
        super_time_start = time.time ()
        
        super_result = [current_task(bias) for bias in small_window]             
        
        super_time_end = time.time ()
        print "Spent %.6f seconds on process %d. \n " % (super_time_end - super_time_start, n)
        return super_result
    def super_task_assemble(results):
        final_results = np.array([])
        for i in range(len(results)):
            final_results = np.append(final_results, results[i])
        
        return np.array( final_results )
    parallel_pool = mp.Pool(processes=4) #automatically uses all cores 
     
    results = parallel_pool.map(super_current_task, [[biaswindow[ (101*n):(101*(n+1))],n] for n in [0, 1, 2, 3]])
    #print results
    results = super_task_assemble(results)
    #print results
    current = results
else:
    current = np.array(experimental_bias/experimental_bias * 1e-10)
###  
scale, error, scaleerror = calculate_error( experimental_bias, current, experimental_current )
print "Error is %2.3e, scale factor %.3e, scale error %.3e" % (error, scale, scaleerror)

minimum = 1.2 * np.min(current)
maximum = 1.2 * np.max(current)

fig = plt.figure(figsize=(12, 10), dpi=1080)
ax = fig.add_subplot(111)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
fig.subplots_adjust(left=0.30)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')

print "Experiment %.3e vs. Theory %.3e" % ( current.max(), experimental_current.max())

delta = lambda VV: np.sqrt( (alpha * VV )**2 + (2. * tau )**2)
epsilon_1 = lambda VV: levels_ni - 0.5 * delta(VV)
epsilon_2 = lambda VV: levels_ni + 0.5 * delta(VV)


analytic_current_0 = lambda VV: realscale * gamma * ( 2. * tau)**2 / (delta(VV)**2 + gamma**2)
analytic_current_1 = lambda VV: np.arctan( (0.5 * VV - epsilon_1(VV))/(gamma/2.0))
analytic_current_2 = lambda VV: np.arctan( (0.5 * VV + epsilon_1(VV))/(gamma/2.0))
analytic_current_3 = lambda VV: np.arctan( (0.5 * VV - epsilon_2(VV))/(gamma/2.0))
analytic_current_4 = lambda VV: np.arctan( (0.5 * VV + epsilon_2(VV))/(gamma/2.0))
analytic_current_5 = lambda VV: gamma/2.0/delta(VV) * np.log( ((0.5*VV-epsilon_1(VV))**2 + (gamma/2.0)**2)/((0.5*VV+epsilon_1(VV))**2 + (gamma/2.0)**2))
analytic_current_6 = lambda VV: gamma/2.0/delta(VV) * np.log( ((0.5*VV-epsilon_2(VV))**2 + (gamma/2.0)**2)/((0.5*VV+epsilon_2(VV))**2 + (gamma/2.0)**2))

analytic_current = lambda VV: analytic_current_0(VV) * ( analytic_current_1(VV)+ analytic_current_2(VV)+ analytic_current_3(VV)+ analytic_current_4(VV)+ analytic_current_5(VV)+ analytic_current_6(VV))

perrin_current = analytic_current( biaswindow)

print "Analytic max %.3e" % perrin_current.max()
print "Theoretical max %.3e" % current.max()
print "Experimental max %.3e" % experimental_current.max()


perrin_scale = experimental_current.max() / perrin_current.max()
perrin_current *= perrin_scale

new_scale = experimental_current.max() / current.max()

current *= new_scale


plt.plot(biaswindow, perrin_current, 'm-', label='non-interacting theoretical') 
plt.plot(biaswindow, current, 'g-', label='interaction theoretical') 
#plt.plot(experimental_bias, original_current, 'b--', label='experimental')   
plt.plot(experimental_bias, experimental_current, 'r--', label='anti-symmetrised experimental')   
#plt.legend()

ax.text( 0.05, 0.40 * experimental_current.min(), "$\\tau=%.3f$" % tau , fontsize=30 )
ax.text( 0.05, 0.52 * experimental_current.min(), "$\\gamma=%.3f$" % gamma , fontsize=30 )
ax.text( 0.05, 0.64 * experimental_current.min(), "$\\alpha=%.3f$" % alpha , fontsize=30 )
ax.text( 0.05, 0.76 * experimental_current.min(), "$\\epsilon_0=%.3f$" % levels , fontsize=30 )
ax.text( 0.05, 0.88 * experimental_current.min(), "$U=%.3f$" % capacitive , fontsize=30 )


plt.legend(bbox_to_anchor=(0., 1.04, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)


title = "Current versus bias"
xlabel = "Bias $V_b$"
ylabel = "$I(V_b)$"
 
#plt.ylim([minimum, maximum])
plt.xlabel(xlabel, fontsize=30)
plt.ylabel(ylabel, fontsize=30)
 
plt.xticks(np.array([-0.25, 0.00, 0.25])) 
plt.yticks(np.linspace( experimental_current.min(), experimental_current.max(), 5))

minorLocator1 = AutoMinorLocator(5)
minorLocator2 = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocator1) 
ax.yaxis.set_minor_locator(minorLocator2) 

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))

plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=20)
plt.tick_params(which='minor', length=10)

plt.title("ni-model %.3e, i-model %.3e" % (1./perrin_scale, 1./new_scale))
#plt.title( "Pts [%s], $\\alpha=%.3f$, $\\tau=%.3f$, $\\Gamma=%.3f$, $\\epsilon_0=%.3f$, $V=%.3f$, $\\beta=%.3f$, $U=%.3f$" % (title,
#    alpha, tau, gamma, levels, bias, beta, capacitive), fontsize=15)     
#plt.legend()


global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)

if plotting_mode == 2 or plotting_mode == 3:
    print "Can't save svg"
    sys.exit(0)
    plt.savefig('perrin_two_site.svg')
else:    
    #plt.show()
    plt.savefig("fit_visualise.png")
