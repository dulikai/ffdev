#! /usr/bin/env python

# force field model: spring case
#
import numpy as np
import math
from bStatus import bStatus


class ffwall:
    """
    a very simple case of particle interaction
    in which,
    spring model: V = 1/2 * epsilon * (r - sigma)^2
    """
    def __init__(self, status):
        self.status = status
        self.params = {'sigma': 2.345, 'epsilon': 1.661}
        self.plane = {}
        return

    def set_plane(self, pA, pB, pC):
        """
        http://keisan.casio.com/exec/system/1223596129
        http://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FEquationofaPlane3Points
        http://www.math.ucla.edu/~ronmiech/Calculus_Problems/32A/chap12/section7/820d37/820_37.html
        known: A, B, C: three points
        ax+by+cz+d=0
        v_n = (a,b,c) = AB \cross AC
        d = -v_n \cdot A
        """
        AB = pB - pA
        AC = pC - pA
        vn = np.cross(AB, AC) # (a,b,c)
        d = - np.dot(vn, A)
        return 
        
    def get_dist(self):
        """
        http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane
        http://mathworld.wolfram.com/Point-PlaneDistance.html
        """
        r = np.linalg.norm(P[0:2])
        x = a * x0 + b * y0 + c * z0 + d
        
        return x / r
        
        
    def energy(self, r):
        """
        V = 1/2 * epsilon * (r - sigma)^2
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
        grad = epsilon * s
        return grad
	
    def interaction(self):
        """
        Write a function that computes the Lennard-Jones potential
        V = 1/2 * epsilon * (r - sigma)^2
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
        
        