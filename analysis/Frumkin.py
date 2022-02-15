import numpy as np
from mkutils import save_to_file, create_fig

# +
# Frumkin model

class FrumkinIsotherm(object):
    '''
    Frumkin Isotherm acc. to Llamas 2018
    https://pubs.acs.org/doi/10.1021/acs.langmuir.8b00358
    '''
    def __init__(self, T=293.15, omega_0=0, b=0, alpha_f=0, epsilon=0, zero_surf=0):
        '''
        Temperature: float, optional
        omega_0: float
            surface area per mol. Unit is m^2/mol.
        b: float
            factor for concentration conversion. Unit is l/mol
        alpha_f: float
            Frumkin parameter, no unit
        epsilon: float
            compressibility of interface (lateral). Unit is m/mN/
        zero_surf: float, optional 
            Interfacial tension without added surfactant in mN/m.
        
        '''
        self.T = T
        self.omega_0 = omega_0
        self.omega_heuristic = omega_0
        self.b = b
        self.alpha_f = alpha_f
        self.epsilon = epsilon
        self.zero_surf = zero_surf

    def concentration(self, gamma):
        '''
        c = omega * gamma / [b*(1 - omega*gamma)] 
            * exp(-2 * alpha_f * omega * gamma)
            
        where gamma is transformed to mol/m2 and omega is iterated
        from:
        
        omega = omega_0 * (1 - epsilon * pi)
        
        for caluclation of pi, see self.get_pi()
        gamma: float 
            molecules per square nanometer (#/nm^2)
            
        Returns: concentration in mol/l
        '''
        # result in mol/l
        if not isinstance(gamma, (list, np.ndarray)):
            gamma = [gamma]
        conc = []
        for gamma_i in gamma:
            omega = self.iterate_omega(gamma_i)
            cap_omega = omega * self.transform_gamma(gamma_i)
            this_conc = (cap_omega / (self.b * (1.-cap_omega)))
            this_conc *= np.exp(-2.*self.alpha_f*cap_omega) 
            conc.append(this_conc)
        return conc
        

    def tension(self, gamma):
        '''
        tens = tens_0  
             + {
             [-R * T / omega] 
             * [ln(1 - omega * gamma) + alpha_f (omega*gamma)^2]
             }
            
        The result in curly brackets is equal to pi.
        where gamma is transformed to mol/m2 and omega is iterated
        from:
        
        omega = omega_0 * (1 - epsilon * pi)
        
        gamma: float 
            molecules per square nanometer (#/nm^2)
            
        Returns: concentration in mol/litre
        '''
        if not isinstance(gamma, (list, np.ndarray)):
            gamma = [gamma]
        tens = []
        for gamma_i in gamma:
            omega = self.iterate_omega(gamma_i)
            pi = self.get_pi(gamma_i, omega)
            tens.append(self.zero_surf - pi)
        return tens

    def iterate_omega(self, gamma):
        '''
        Iterate omega = omega_0 * (1 - epsilon*pi(omega)).
        Start with a guess/previous result (self.omega_heuristic).
        '''
        omega = self.omega_heuristic
        d_omega = 1
        max_omega = self.omega_0 * (1.0 - self.epsilon * self.zero_surf)
        i = 0

        while d_omega > 1e-6 and i < 1000:
            pi = self.get_pi(gamma, omega)
            new_omega = self.omega_0 * (1.0 - (self.epsilon * pi))
            d_omega = abs(omega - new_omega)
            # print(cap_omega, pi, omega, new_omega)
            omega = new_omega
            i += 1
        if not np.isnan(omega):
            self.omega_heuristic = omega
        return omega
        
        

    def get_pi(self, gamma, omega):
        '''    
        pi = [-R * T / omega] 
             * [ln(1 - omega * gamma) + alpha_f (omega*gamma)^2]
        '''    
        RT = 8.3145 * self.T
        #cap_omega = min(omega * self.transform_gamma(gamma), 1.0)
        cap_omega = omega * self.transform_gamma(gamma)
        pi = (
            -RT
            / omega
            * (np.log(1.0 - cap_omega) + (self.alpha_f * cap_omega * cap_omega))
            * 1000.
        )
        return pi

    @staticmethod
    def transform_gamma(gamma):
        '''
        from 1/nm^2 to mol/m^2
        constant: 10^18/(6.022*10^23) nm2/m2 * mol/#
        '''
        return gamma / (6.022e5)

    def get_gamma_zero_surftens(self):
        '''
        omega in m2/mol, transform with nm2/m2 * mol/#
        to nm2/#, then inverse
        '''
        omega = self.get_omega_zero_surftens()
        return 1.0 / (omega * self.transform_gamma(1))
    
    def get_omega_zero_surftens(self):
        '''
        return omega when pi = zero_surftens and therefore 
        the interfacial tension is zero.
        '''
        return self.omega_0 * (1.0 - self.epsilon * self.zero_surf)