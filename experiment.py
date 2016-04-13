import numpy as np
import scipy.interpolate as si
def read_experiment_data( filename ):
    
    file_handler = open( filename, "r" );
    
    data = np.genfromtxt(file_handler, skip_header=9, dtype=None, usecols=range(0,3)); #excluding the symtype col
    exp_bias = data[:,0]
    exp_current = data[:,1] 
    
    truebias = np.linspace(np.min(exp_bias), np.max(exp_bias), int(exp_bias.shape[0]/2))
    
    truecurrent = si.griddata(exp_bias, exp_current, truebias, method='nearest')
    file_handler.close()
    
    return truebias, truecurrent
def read_experiment( filename ):
    exp_bias, exp_current = read_experiment_data(filename)
    
    dI = np.max(exp_current) - np.min(exp_current)
    dV = np.max(exp_bias) - np.min(exp_bias)
    
    exp_background = lambda V: dI/dV * V
    
    exp_current -= exp_background(exp_bias)
    
    return exp_bias, exp_current 
