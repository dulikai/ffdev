#! /usr/bin/env python

# force field model: linear case

# lj for trans. and 
# spin for rot.
#
import numpy as np
import math
from bStatus import bStatus


class fflinear:
    """
    a very simple case of particle interaction
    in which,
    LJ for com, self spin for rot.
    """
    def __init__(self, status):
        self.status = status
        self.params = {'sigma': 2.345, 'epsilon': 1.661}
        return

    def energy(self, r):
        """
        V_{rep-disp} = 4 \epsilon \cdot [V_{rep} - V_{disp}]
        x = sigma/r; V_{rep} = x^12; V_{disp} = x^6
        """
        sigma = self.params['sigma']
        epsilon = self.params['epsilon']
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        pot = 4.0 * epsilon * (s12 - s6)
        return pot
        
    def gradient(self, r):
        """
        V' = 4 \epsilon [-12/r x^12 - (-6/r) x^6]
        """
        sigma = self.params['sigma']
        epsilon = self.params['epsilon']
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        grad = 4.0 * epsilon * ((-12.0/r) * s12 - (-6/r) * s6)
        grad = 0.5 * (r - 5.0)
        return grad
	
    def trans_eandg(self):
        """
        Write a function that computes the Lennard-Jones potential
        V = 4 \epsilon [(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6]
        dist = self.get_dist_mat()
        """
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
        # for ibody in sites:
            # print "na force", ibody.force
        return
        
    def rot_eandg(self):
        """ rot. mod """
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        for body in sites:
            body.torque = np.random.random(3)
        self.status.potential_energy += np.random.random(1)
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
        
    def nd_gradient(self):
        """
        numerical force
        """
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
            # print "nd force", ibody.force
        return
                
    def eandg_nd(self):
        self.status.zeros()
        ene = self.compute_energy()
        self.status.potential_energy = ene
        self.nd_gradient()
        self.rot_eandg()
        return
        
    def eandg(self):
        """
        energy & gradient
        """
        self.status.zeros()        
        self.trans_eandg()
        self.rot_eandg()
        # print self.status.potential_energy
        return
        
        
if __name__ == "__main__":
    status = bStatus()
    ff = fflinear(status)
    ff.eandg()
    ff.nd_gradient()
        
        