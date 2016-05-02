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

plotting_mode = 0

alpha = 0.74
tau = -0.0241
gamma = 0.0102
levels = -0.25
bias = 0.25*4


#levels set to zero-
levels = -0.05

hamiltonian = np.zeros((2,2))

hamiltonian[0][0] = levels + 0.5 * alpha * bias
hamiltonian[1][1] = levels - 0.5 * alpha * bias

tunnel = np.zeros((2,2))
tunnel[0][1] = -tau
tunnel[1][0] = -tau

#change these numbers based on visual inspection
epsilon_left = -.5
epsilon_right = 1.
epsilon_res = 400 


gamma_left = np.zeros((2,2))
gamma_left[0][0] = gamma

gamma_right = np.zeros((2,2))
gamma_right[1][1] = gamma

beta = 250.0

epsilon = np.linspace(epsilon_left, epsilon_right, epsilon_res);

data_capacitive = []
data_levels = []
data_tau = []

for capacitive in np.arange(0, 0.5, 0.01):
    
    interaction = np.zeros((2,2))
    interaction[0][1] = capacitive
    interaction[1][0] = capacitive
    
    calculation = igfwl(
        hamiltonian, 
        tunnel,
        interaction, 
        gamma_left,
        gamma_right, 
        beta
    )
    
    transmission = calculation.full_transmission(epsilon)
    
    ##simple calculation for non-interacting
    el = hamiltonian[0, 0]
    er = hamiltonian[1, 1]
    
    delta = np.sqrt( (el-er)**2 + (2*tau)**2)
    
    left_eig = 0.5*(el+er) - 0.5 * delta
    right_eig = 0.5*(el+er) + 0.5 * delta
    
    def peakfinder(tg, e): 
        """Looks for the transmission peaks given the gradient and energy, assuming there are but two peaks."""
        peaks = e[tg>np.max(tg)*0.5]
        valleys =  e[tg<np.min(tg)*0.5]
        
        left = []
        right = []
        
        peak_average = np.average(peaks)
        for p in peaks:
            if p < peak_average:
                left.append(p)
            else:
                right.append(p)
                
        valley_average = np.average(valleys)
        for v in valleys:
            if v < valley_average:
                left.append(v)
            else:
                right.append(v)
        
        
        
        left_mu = np.average(left)
        left_sigma = np.std(left)

        right_mu = np.average(right)
        right_sigma = np.std(right)

        return np.array([[left_mu, left_sigma],[right_mu, right_sigma]])        
        
    grad = np.gradient(transmission)
    
    located = peakfinder(grad, epsilon)
    
    lmin = located[0,0]
    lmax = located[1,0]
    
    effective_levels = (lmin + lmax)/2
    delta = lmax - lmin
    
    effective_tau = -np.sqrt(-0.25*(alpha**2 * bias**2-delta**2))
    
    print capacitive, effective_levels, effective_tau
    
    data_capacitive.append(capacitive)
    data_levels.append(effective_levels)
    data_tau.append(effective_tau)
    
    
    
plt.figure(figsize=(10, 10), dpi=1080)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
 
 
title = "Effective parameters due to interaction. $\\beta=%.3f$, $V=%.3f$, $\\epsilon_0=%.3f$, $\\Gamma=%.3f$, $\\tau=%.3f$,$\\alpha=%.3f$, " %(beta,
    bias, levels, gamma, tau, alpha)
xlabel = "Capacitive Interaction $U$"
ylabel = "Strength (eV)"
plt.rc('font', family='serif')

data_capacitive = np.array(data_capacitive)
data_levels = np.array(data_levels)
data_tau = np.array(data_tau)

taurico = -(data_tau.max() - data_tau.min())/0.45
print "dtau = %.3f for V=%.3f" % (taurico, bias)

plt.plot(data_capacitive, data_levels, 'go', label='effective $\\epsilon_0$')    
plt.plot(data_capacitive, data_tau, 'ro', label='effective $\\tau$')  


plt.xlabel(xlabel, fontsize=15)
plt.ylabel(ylabel, fontsize=15)

plt.title( "%s" % (title), fontsize=15)     
plt.legend(loc='upper left')

plt.show() 