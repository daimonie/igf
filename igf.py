import numpy as np
from numba import jit

class igfwl(object): 
    """The igfwl class will compute the transport function using the many-body interacting Green's function in the wide-band limit."""
    def __init__(self, 
        param_epsilon, 
        param_tau,
        param_u, 
        param_gamma_left,
        param_gamma_right,
        param_beta): 
            """The init function takes the single-particle energy matrix (diagonal), the inter-orbital hopping(tunneling), the interaction energy, and the coupling to the leads (Gamma). Finally, it holds the inverse temperature (initial) and the fermi energies of the leads."""
            self.epsilon     = param_epsilon
            self.tau         = param_tau
            self.u           = param_u
            self.gamma_left  = param_gamma_left
            self.gamma_right = param_gamma_right
            
            self.sigma_retarded = 1j * (self.gamma_left + self.gamma_right)
            self.sigma_advanced = - self.sigma_retarded;
            
            self.dim    = len(self.u)
            self.rho = np.zeros((2**self.dim))
            
            self.beta = param_beta
            
    # jit decorator tells Numba to compile this function.
    # The argument types will be inferred by Numba when function is called.
    @jit     
    def singleparticlebackground(self, background):
        """This gives us the single-particle Green's function with some background."""
        
        mu_background = np.diag([self.epsilon[i][i] + np.dot( self.u[i], background) for i in range(0,self.dim)])
        
        single_retarded = lambda energy: np.linalg.inv( np.eye( self.dim) * energy - mu_background - self.tau - self.sigma_retarded)
        single_advanced = lambda energy: np.linalg.inv( np.eye( self.dim) * energy - mu_background - self.tau - self.sigma_advanced)
        
        return single_retarded, single_advanced    
    # jit decorator tells Numba to compile this function.
    # The argument types will be inferred by Numba when function is called.
    @jit     
    def generate_superset(self, number):
        """ This function returns a list of the integers that are in the superset of k. The reason this is a seperate function is because I think this algorithm can be a lot smoother/faster."""
        
        superset = []
        for i in range(0, 2**(self.dim)):
            #Note: number == i is excluded by the P_kk P_ll rule anyway.
            #if  i != number and (number & i)==number:
            if  (number & i)==number:
                superset.append(i)
        return superset    
    # jit decorator tells Numba to compile this function.
    # The argument types will be inferred by Numba when function is called.
    @jit     
    def number(self, ket):
        """Returns the number of the ket"""
        
        final = 0.0
        q = 0
        for i in ket:
            if i != 0:
                final += 2**q
            q += 1
            
        return final    
    # jit decorator tells Numba to compile this function.
    # The argument types will be inferred by Numba when function is called.
    @jit     
    def ket(self, number):
        """Turns an integer number into a ket."""
        
        final_ket = np.array( [(2**i==number&2**i)*1.0 for i in range(0,self.dim)] )
        return final_ket     
    # jit decorator tells Numba to compile this function.
    # The argument types will be inferred by Numba when function is called.
    @jit     
    def distribution(self):
        """Sets the boltzmann distribution function and generates the density-matrix."""
        
        
        energy_vector = []
        superset = self.generate_superset(0) 
        
        for i in superset:
            state           = self.ket(i)
            
            norm_squared    = np.dot(state.T, state)
            
            if norm_squared > 0: #zero is appended at the end
                energy          = np.dot(state.T, np.dot( self.epsilon, state)) / norm_squared
                interaction     = np.dot(state.T, np.dot( self.u, state)) / norm_squared
                 
                energy_vector.append( energy + interaction )
                
        energy_vector.insert(0, 0.0)
        
        probability = np.exp( np.multiply(-self.beta, energy_vector))
        probability /= probability.sum()
        
        return probability    
    # jit decorator tells Numba to compile this function.
    # The argument types will be inferred by Numba when function is called.
    @jit     
    def transport_channel(self, k, epsilon):
        """Returns the transmission function for the many body state k."""
        advanced_gf = np.zeros((self.dim, self.dim))
        retarded_gf = np.zeros((self.dim, self.dim))
        
        transport_k = 0 * epsilon
        
        chances = self.distribution()
        
        for i in self.generate_superset(k):
                state = self.ket( i ) 
                chance = chances[i]  
                ret_gf, ad_gf = self.singleparticlebackground( state ) 
                   
            
                transport_k_ij = [np.trace(chance**2 *np.dot(self.gamma_left, ( np.dot(
                    ret_gf(ee),  np.dot(self.gamma_right, ad_gf(ee)))))) for ee in epsilon]
                transport_k += np.real(transport_k_ij)
        return transport_k    
    # jit decorator tells Numba to compile this function.
    # The argument types will be inferred by Numba when function is called.
    @jit     
    def spectral_channel(self, k, epsilon):
        """Returns the transmission function for the many body state k."""
        advanced_gf = np.zeros((self.dim, self.dim))
        retarded_gf = np.zeros((self.dim, self.dim))
        
        transport_k = 0 * epsilon
        
        chances = self.distribution()
        
        for i in self.generate_superset(k):
                state = self.ket( i ) 
                
                ret_gf, ad_gf = self.singleparticlebackground( state ) 
                
                chance = chances[i] 
            
                transport_k_ij = [np.trace(chance**2 * 1.0/np.pi * np.imag(ret_gf(ee))) for ee in epsilon]
                transport_k += np.real(transport_k_ij)
        return transport_k