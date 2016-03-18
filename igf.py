import numpy as np
import scipy as sc
import scipy.misc as scmisc
import scipy.special as scspecial

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
            
            self.sigma_retarded = 1j * (self.gamma_left + self.gamma_right) / 2.0
            self.sigma_advanced = - self.sigma_retarded;
            
            self.dim    = len(self.u)
            self.rho = np.zeros((2**self.dim))
            
            self.beta = param_beta
    
            self.cutoff_chance = 0.0001
            
    def singleparticlebackground(self, background):
        """This gives us the single-particle Green's function with some background."""
        
        mu_background = np.diag([self.epsilon[i][i] + np.dot( self.u[i], background) for i in range(0,self.dim)])
        
        single_retarded = lambda energy: np.linalg.inv( np.eye( self.dim) * energy - mu_background - self.tau - self.sigma_retarded)
        single_advanced = lambda energy: np.linalg.inv( np.eye( self.dim) * energy - mu_background - self.tau - self.sigma_advanced)
        
        return single_retarded, single_advanced
    def generate_superset(self, number):
        """ This function returns a list of the integers that are in the superset of k. The reason this is a seperate function is because I think this algorithm can be a lot smoother/faster."""
        
        superset = []
        for i in range(0, 2**(self.dim)):
            #Note: number == i is excluded by the P_kk P_ll rule anyway.
            #if  i != number and (number & i)==number:
            if  (number & i)==number:
                superset.append(i)
        return superset
    def number(self, ket):
        """Returns the number of the ket"""
        
        final = 0.0
        q = 0
        for i in ket:
            if i != 0:
                final += 2**q
            q += 1 
        return final
    def ket(self, number):
        """Turns an integer number into a ket."""
        
        final_ket = np.array( [(2**i==number&2**i)*1.0 for i in range(0,self.dim)] )
        return final_ket 
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
    
    def transport_channel_ij(self, i, chances, epsilon):
        state = self.ket( i ) 
        chance = chances[i]  
        ret_gf, ad_gf = self.singleparticlebackground( state ) 
         
        transport_k_ij = [np.trace(chance**2 *np.dot(self.gamma_left, ( np.dot(
            ret_gf(ee),  np.dot(self.gamma_right, ad_gf(ee)))))) for ee in epsilon]
        return transport_k_ij
    def transport_channel(self, k, epsilon):
        """Returns the transmission function for the many body state k."""
        advanced_gf = np.zeros((self.dim, self.dim))
        retarded_gf = np.zeros((self.dim, self.dim))
        
        transport_k = 0 * epsilon
        
        chances = self.distribution()
        for i in self.generate_superset(k):
                if chances[i] > self.cutoff_chance:
                    transport_k += np.real( self.transport_channel_ij(i, chances, epsilon))
        return transport_k
    def spectral_channel(self, k, epsilon):
        """Returns the transmission function for the many body state k."""
        advanced_gf = np.zeros((self.dim, self.dim))
        retarded_gf = np.zeros((self.dim, self.dim))
        
        transport_k = 0 * epsilon
        
        chances = self.distribution()
        
        superset = self.generate_superset(k)
        
        for i in superset:
                state = self.ket( i ) 
                
                ret_gf, ad_gf = self.singleparticlebackground( state ) 
                
                chance = chances[i] 
            
                transport_k_ij = [np.trace(chance**2 * 1.0/np.pi * np.imag(ret_gf(ee))) for ee in epsilon]
                transport_k += np.real(transport_k_ij) / len(superset)
        return transport_k
#################
class igfwl_vibrational(igfwl):
    def __init__(self, 
        param_epsilon, 
        param_tau,
        param_u, 
        param_gamma_left,
        param_gamma_right,
        param_beta,
        phonon_energy,
        phonon_electron_coupling ): 
        #call parent init
        super(igfwl_vibrational, self).__init__( 
            param_epsilon, 
            param_tau,
            param_u, 
            param_gamma_left,
            param_gamma_right,
            param_beta)
        #own
        self.pe = phonon_energy
        self.pec = phonon_electron_coupling
        
    def overlap(self, n, m):
        n = int(n)
        m = int(m)
        
        l = self.pec
        
        if n >= 0 and m >= 0:
            smaller =  np.min( [n, m] )
            
            fc = np.sqrt(scspecial.factorial(n) * scspecial.factorial(m))
            if np.isinf(fc):
                raise Exception("[igfwl_vibrational::overlap] Numbers are too high. Can't cope.")
            
            fc *= np.exp(-0.5 * l**2)
            
            fc *= np.sum([ (-l)**(n-k) * l**(m-k) / scspecial.factorial(k) / scspecial.factorial(m-k) / scspecial.factorial(n-k)   for k in range(smaller+1)])
            
            fc = fc
            return fc
        return 0.0
    def sum_factor(self, n, m, x, y): 
        
        l = self.pec
        
        fc = self.overlap(n+y-x, m-y)
        sc = (l * self.pe)**x  
        sc *= scmisc.comb(x, y)
        
        sc *= fc 
        
        return sc
    def chance_phonons(self, m, n):
        p_m = np.exp(-self.beta * self.pe * (m+n)/2.0) * (1 - np.exp(-self.beta * self.pe))
        
        return p_m
    def transport_channel_vibrational(self, k, epsilon, significant_terms):
        """Returns the transmission function for the many body state k.""" 
        transport_k = []
        
        chances = self.distribution()
        
        for energy in epsilon:
            retA_vib = 0
            advA_vib = 0
            for i in self.generate_superset(k):
                if chances[i] > self.cutoff_chance:
                    retA, advA = self.singleparticlebackground(i)
                    for significant_term in significant_terms:   
                        retA_vib += retA(energy)**significant_term[2] * significant_term[4]
                        advA_vib += advA(energy)**significant_term[2] * significant_term[4]
            transport_at_energy = np.abs( np.real(np.trace(chances[i]**2 *np.dot(self.gamma_left, ( np.dot(retA_vib,  np.dot(self.gamma_right, advA_vib)))))) )
            transport_k.append(transport_at_energy)
            print "%d\t%2.3f\t%2.3e" % (i,energy,transport_at_energy)
        return np.array(transport_k )