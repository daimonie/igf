import sys as sys 
import numpy as np
from paralleltask import *
import scipy.interpolate as si
from scipy.optimize import minimize
from scipy.constants import physical_constants as pc
from igf import *
#Command line arguments.
import argparse as argparse  
import time
from experiment import *
##matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm
from scipy.constants import physical_constants as pc
###
global_time_start = time.time()
def kitty():
    print "Mew"
    sys.exit(0)
###

parser  = argparse.ArgumentParser(prog="roughfit visualise",
  description = "Plots experimental current, calculated current and finally the analytical non-interacting current.")  
  
parser.add_argument(
    '-s',
    '--sep',
    help='Separation data file',
    action='store',
    type = int,
    default = 638
) 
parser.add_argument(
    '-a',
    '--alpha', 
    action='store',
    type = float,
    default = .75
)    
parser.add_argument(
    '-t',
    '--tau', 
    action='store',
    type = float,
    default = .01
)    
parser.add_argument(
    '-g',
    '--gamma', 
    action='store',
    type = float,
    default = .001
)    
parser.add_argument(
    '-e',
    '--epsilon', 
    action='store',
    type = float,
    default = -.20
)    
parser.add_argument(
    '-u',
    '--capacitive', 
    action='store',
    type = float,
    default = .2
)     
parser.add_argument(
    '-z',
    '--zeta', 
    action='store',
    type = float,
    default = 1.00
)    
parser.add_argument(
    '-k',
    '--ksi', 
    action='store',
    type = float,
    default = 1.00
)  
parser.add_argument(
    '-c',
    '--cores', 
    action='store',
    type = int,
    default = 1
)    
parser.add_argument(
    '-m',
    '--mode', 
    action='store',
    type = int,
    default = 1
)    

args    = parser.parse_args()  
sep = args.sep

tau         = args.tau
gamma       = args.gamma
levels      = args.epsilon - 1e-4
alpha       = args.alpha
capacitive  = args.capacitive 
cores       = args.cores
zeta        = args.zeta
ksi         = args.ksi
mode        = args.mode
beta        = 250.00
if mode < 0 or mode > 3:
    raise Exception("Improper mode.")
###

separation_array = range(638, 670)        
        
data_bias = np.zeros(( len(separation_array)-1, 404))
data_current = np.zeros(( len(separation_array)-1, 404))

i = 0          
for sep in separation_array:
    if sep != 645:
        file = "exp_data/IV130328_7_%d.dat" % sep
        #print "Reading [%s]" % file
        bias, current = read_experiment(file)
     
        filter = np.ones(15)/15.0
        current = np.convolve(current, filter, mode='same')
         
        bias_array = bias 
                
        data_bias[i]    = bias
        data_current[i] = current
        
        i += 1
###
bias         = data_bias.mean(axis=0)
experimental = data_current.mean(axis=0) / 1e-9
### 
def calculate_formula(arguments): 
    epsilon_res = 250
    bias         = arguments[0]
    alpha       = arguments[1]
    tau         = arguments[2]
    gamma       = arguments[3]
    capacitive  = arguments[4]
    beta        = arguments[5]
    levels      = arguments[6]
    levels_ni   = levels
    ksi         = arguments[7]
    realscale   = pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]
    
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

    current = analytic_current( bias)
    return [bias, current]

def calculate_spinfull(arguments): 
    epsilon_res = 250
    bias        = arguments[0]
    alpha       = arguments[1]
    tau         = arguments[2]
    gamma       = arguments[3]
    capacitive  = arguments[4]
    beta        = arguments[5]
    levels      = arguments[6]
    ksi         = arguments[7]
    zeta        = arguments[8]
    levels_ni   = levels
    realscale   = pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]
    
    hamiltonian = np.zeros((4,4))

    hamiltonian[0][0] = levels + 0.5 * alpha * bias
    hamiltonian[1][1] = levels + 0.5 * alpha * bias
    hamiltonian[2][2] = levels - 0.5 * alpha * bias
    hamiltonian[3][3] = levels - 0.5 * alpha * bias

    tunnel = np.zeros((4,4))
    tunnel[0][2] = tau  
    tunnel[1][3] = tau  

    #right-left
    tunnel[2][0] = tau  
    tunnel[3][1] = tau   

    interaction = np.zeros((4,4))

    interaction[0][1] = ksi * capacitive
    interaction[1][0] = ksi * capacitive

    interaction[2][3] = ksi * capacitive
    interaction[3][2] = ksi * capacitive

    interaction[0][3] = zeta*capacitive
    interaction[3][0] = zeta*capacitive

    interaction[0][2] = zeta*capacitive
    interaction[2][0] = zeta*capacitive

    interaction[1][3] = zeta*capacitive
    interaction[3][1] = zeta*capacitive

    interaction[1][2] = zeta*capacitive
    interaction[2][1] = zeta*capacitive

    ##print interaction
    #sys.exit(0)
    gamma_left = np.zeros((4,4))
    gamma_left[0][0] = gamma
    gamma_left[1][1] = gamma

    gamma_right = np.zeros((4,4))
    gamma_right[2][2] = gamma
    gamma_right[3][3] = gamma
    
    calculation = igfwl(
        hamiltonian, 
        tunnel,
        interaction, 
        gamma_left,
        gamma_right, 
        beta
    ) 
    
    epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
   
    spinfull_transmission = calculation.full_transmission(epsilon)  
    spinfull_current = realscale*np.trapz(spinfull_transmission, epsilon)

    return [bias, spinfull_current]
def calculate_spinless(arguments):  
    epsilon_res = 250
    bias        = arguments[0]
    alpha       = arguments[1]
    tau         = arguments[2]
    gamma       = arguments[3]
    capacitive  = arguments[4]
    beta        = arguments[5]
    levels      = arguments[6]
    ksi         = arguments[7]
    zeta        = arguments[8]
    levels_ni   = levels
    realscale   = pc["elementary charge"][0] / pc["Planck constant"][0] * pc["electron volt"][0]

    spinless_hamiltonian = np.zeros((2,2))
    spinless_hamiltonian[0][0] = levels + 0.5 * alpha * bias
    spinless_hamiltonian[1][1] = levels - 0.5 * alpha * bias

    spinless_tunnel = np.zeros((2,2))
    spinless_tunnel[0][1] = -tau
    spinless_tunnel[1][0] = -tau

    spinless_interaction = np.zeros((2,2))
    spinless_interaction[0][1] = capacitive
    spinless_interaction[1][0] = capacitive


    spinless_gamma_left = np.zeros((2,2))
    spinless_gamma_left[0][0] = gamma

    spinless_gamma_right = np.zeros((2,2))
    spinless_gamma_right[1][1] = gamma

    spinless_calculation = igfwl(
        spinless_hamiltonian, 
        spinless_tunnel,
        spinless_interaction, 
        spinless_gamma_left,
        spinless_gamma_right, 
        beta
    )

    
    epsilon = np.linspace(-bias/2.0, bias/2.0, epsilon_res);
  
    spinless_transmission = spinless_calculation.full_transmission(epsilon)
    spinless_current = realscale*np.trapz(spinless_transmission, epsilon) 
     
    
    return [bias, spinless_current]


#### calculate current
manager = 0
if mode == 0:
    manager = taskManager( cores, calculate_formula ) 
elif mode == 1:
    manager = taskManager( cores, calculate_spinless ) 
elif mode == 2:
    manager = taskManager( cores, calculate_spinfull ) 

for this_bias in bias:      
#for this_bias in np.linspace(-1.0, 1.0, 100):
    manager.add_params([this_bias, alpha, tau, gamma, capacitive, beta, levels, ksi, zeta]) 
    #print [this_bias, alpha, tau, gamma, capacitive, beta, levels, ksi, zeta]
    #sys.exit(0)
manager.execute()
results = manager.final()
results = np.array(results)

calculated_bias = results[:,0]
calculated_current = results[:,1]/1e-9

#######################################
fig = plt.figure(figsize=(20, 20), dpi=1080)
ax = fig.add_subplot(111)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
fig.subplots_adjust(left=0.30)

title = "Dummy title"
xlabel = ""
ylabel = ""
plt.rc('font', family='serif')


print "Max current %2.3e" % np.max(np.abs(calculated_current))

make_same_scale = np.max(experimental) / np.max(calculated_current)
calculated_current *= make_same_scale

plt.plot(bias, experimental, 'm-', label='Experimental Average') 
if mode == 0: #formula
    plt.plot(calculated_bias, calculated_current, 'r-', label='Non-Interacting two-site model') 
elif mode == 1: #spinless
    plt.plot(calculated_bias, calculated_current, 'g-', label='Interacting spinless model') 
elif mode == 2: #spinfull
    plt.plot(calculated_bias, calculated_current, 'b-', label='Interacting spinfull model') 
#plt.legend()
if mode == 2:
    ax.set_title("Scaled $I(V)$ by $%.2f$, $\\tau=%.3f, \\gamma=%.3f, \\alpha=%.3f, \\epsilon_0=%.3f, U=%.3f, \\xi = 1.0, \\zeta = 1.0$" % (make_same_scale, tau, gamma, alpha, levels, capacitive), fontsize=20)
else:
    ax.set_title("Scaled $I(V)$ by $%.2f$, $\\tau=%.3f, \\gamma=%.3f, \\alpha=%.3f, \\epsilon_0=%.3f, U=%.3f$" % (make_same_scale, tau, gamma, alpha, levels, capacitive), fontsize=20)

plt.legend(bbox_to_anchor=(0., 1.04, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0., fontsize=25)

 
xlabel = "Bias $V_b$ [eV]"
ylabel = "Current $I(V_b)$  [nA] "
 
plt.xlim([-0.25, 0.25])
#plt.xlim([-1., 1.])
plt.xlabel(xlabel, fontsize=25)
plt.ylabel(ylabel, fontsize=25)
 
plt.xticks(np.array(range(11))*0.05-0.25) 

minorLocator1 = AutoMinorLocator(5)
minorLocator2 = AutoMinorLocator(5)
ax.xaxis.set_minor_locator(minorLocator1) 
ax.yaxis.set_minor_locator(minorLocator2) 

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=20)
plt.tick_params(which='minor', length=10)

plt.savefig('fit_average.pdf')
###
global_time_end = time.time ()
print "\n Time spent %.6f seconds. \n " % (global_time_end - global_time_start)