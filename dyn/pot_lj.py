#! python

# Lennard John potential 
import numpy as np

import pgrad as pg

class pot_lj126():
	""" 
	cluster summary
	loop for... 
	"""
    def __init__(self, model):
        return
    def energy(self):
        """
        V_{rep-disp} = 4 \epsilon \cdot [V_{rep} - V_{disp}]
        x = sigma / r
        V_{rep} = x^12
        V_{disp} = x^6
        """
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        pot = 4.0 * epsilon * (s12 - s6)
        return pot
        
    def gradient(self):
        """
        V' = 4 \epsilon [-12/r x^12 - (-6/r) x^6]
        """
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        grad = 4.0 * epsilon * ((-12.0/r) * s12 - (-6/r) * s6)
        return grad
	
    def eandg(self, r, sigma, epsilon):
        """
        energy & gradient
        """
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        ene = 4.0 * epsilon * (s12 - s6)
        grad = 4.0 * epsilon * ((-12.0/r) * s12 - (-6/r) * s6)
        return ene, grad
        
	def calc(self):
		"""
		calculate the total e & g
		"""
		

        rij = rj - ri
        \frac{\partial r_{ij}}{\partial r_k} = 
        (\delta_{jk} - \delta_{ik}) \frac{\vec{\bf r}_{ij}}{r_{ij}}
        r = np.linalg.norm(rij)
        i_deta = dirac_deta(j,k) - dirac_deta(i,k);
        ur = rij / r * i_deta
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        ene = 4.0 * epsilon * (s12 - s6)
        grad = 4.0 * epsilon * ((-12.0/r) * s12 - (-6/r) * s6)
		
    def get_dist_mat(self):
        """
        calculate the dist matrix.
        """
        sites = self.md['site']
        n_site = self.md['n_sites']
        self.dist_mat = [[0.0 for j in xrange(n_site)] for i in xrange(n_site)]
        for i in xrange(n_site):
            for j in xrange(n_site):
                coordi = site[i]['coord']
                coordj = site[j]['coord']
                self.dist_mat[i][j] = np.linalg.norm(coordi-coordj)
        return
    def get_potential_energy(self):
        """
        Write a function that computes the Lennard-Jones potential
        V = 4 \epsilon [(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6]
        """
        lj.get_distances()
        sigma = 0.3345
        epsilon = 0.0661
        dist_mat = self.dist_mat
        n_site = self.cluster['n_site']
        v = 0.0
        for i in xrange(n_site):
            for j in xrange(i+1, n_site):
                r = dist_mat[i][j]
                s = sigma / r
                s6 = s**6; s12 = s6*s6;
                v += 4.0 * epsilon * (s12 - s6)
        return v
        
        
        
        
        