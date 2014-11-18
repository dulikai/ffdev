#! /usr/bin/env python
#
# ----------------------------------------------------
# @ dulikai
# @ 2014.11.18
# @ qibebt
# @ two state model
#
# ----------------------------------------------------
# force field model for test only, so simple
# dulikai
# @2014.10.18
# @ qibebt
# 
# harmonic potential for inter site force
# ____________________________________________________
#

import os
import numpy as np
import math
from bStatus import bStatus


class ffTest:
    """
    a very simple case of particle interaction
    in which,
    harmonic potential
    """
    def __init__(self, status):
        self.status = status
        self.params = {'sigma1': 2.345, 'epsilon1': 1.661,
                       'sigma2': 3.500, 'epsilon2': 2.500}
        return

    # v11 = 0.5 * k_11 * (r - r_11)
    # v22 = 0.5 * k_22 * (r - r_22)
    # v12 = v21 = 1.0
    def hamiltonian(self):
        """ two state model """
        r_11 = self.params['sigma1']
        k_11 = self.params['epsilon1']
        r_22 = self.params['sigma2']
        k_22 = self.params['epsilon2']
        # v_11
        s = r - r_11
        v_11 = 0.5 * k_11 * s * s
        # v_22 
        s = r - r_22
        v_22 = 0.5 * k_22 * s * s
        # v_12, v_21
        v_12 = 1.0
        v_21 = v_12
        #
        H = np.zeros((2,2))
        H[0][0] = v_11; H[0][1] = v_12
        H[1][0] = v_21; H[1][1] = v_22
        
        return H
        
    def gradient(self):
        """
        V' = epsilon * (r - sigma)
        """
        sigma = self.params['sigma']
        epsilon = self.params['epsilon']
        s = r - sigma
        g = epsilon * (r - sigma)
        return g
    
    def energy(self, r):
        """
        V = 0.5 * epsilon * (r-sigma)^2
        """
        sigma = self.params['sigma']
        epsilon = self.params['epsilon']
        s = r - sigma
        pot = 0.5 * epsilon * s * s
        return pot
        
    def gradient(self, r):
        """
        V' = epsilon * (r - sigma)
        """
        sigma = self.params['sigma']
        epsilon = self.params['epsilon']
        s = r - sigma
        g = epsilon * (r - sigma)
        return g
	
    def eandg(self):
        """
        Write a function that computes the Lennard-Jones potential
        V = 0.5 * epsilon * (r-sigma)^2
        dist = self.get_dist_mat()
        """
        self.status.zeros()  
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        # print dist
        for i in xrange(n_site):
            # force act on site i. f_r dr_i
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos -ibody.pos
                r = np.linalg.norm(rij)
                drdi = - rij / r
                drdj = - drdi
                self.status.potential_energy += self.energy(r)    
                ibody.force += - self.gradient(r) * drdi
                jbody.force += - self.gradient(r) * drdj
        for ibody in sites:
            print "na force", ibody.force
        return
        
        
        
    def get_dist_mat(self):
        """
        calculate the distance matrix.
        """
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        dist_mat = [[0.0 for j in xrange(n_site)] for i in xrange(n_site)]
        for i in xrange(n_site):
            for j in xrange(n_site):
                ri = sites[i].pos
                rj = sites[j].pos
                dist_mat[i][j] = np.linalg.norm(ri-rj)
                # print ri, rj
        return dist_mat        
        
    def compute_energy(self):
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        # print dist
        ene = 0.0
        for i in xrange(n_site):
            # force act on site i. f_r dr_i
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos -ibody.pos
                r = np.linalg.norm(rij)
                ene += self.energy(r)    
        return ene
        
    def eandg_nd(self):
        """
        numerical force
        """
        self.status.zeros()
        # energy
        ene = self.compute_energy()
        self.status.potential_energy = ene
        # grad
        delta = 1.0e-6
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        for ibody in sites:
            for j in xrange(3):
                ibody.pos[j] += delta
                v2 = self.compute_energy()
                ibody.pos[j] += -2.0 * delta
                v1 = self.compute_energy()
                ibody.pos[j] += delta
                ibody.force[j] = - (v2-v1) / (2.0 * delta)
            print "nd force", ibody.force
        return
                
     
     
# main program        
if __name__ == "__main__":
    os.chdir("test")
    status = bStatus()
    ff = ffTest(status)
    ff.eandg()
    ff.eandg_nd()
        
        
        
        
        