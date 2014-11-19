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


class ffTwoState:
    """
    a very simple case of particle interaction
    in which,
    harmonic potential
    """
    def __init__(self, status):
        self.status = status
        # sigma epsilon pairs
        self.params = [{'sigma': 2.345, 'epsilon': 1.661}, # spin 1
                       {'sigma': 3.500, 'epsilon': 2.500}]  # spin 2
        return

        
    def energy(self, r, sigma, epsilon):
        """
        V = 0.5 * epsilon * (r-sigma)^2
        """
        s = r - sigma
        pot = 0.5 * epsilon * s * s
        return pot
        
    def gradient(self, r, sigma, epsilon):
        """
        V' = epsilon * (r - sigma)
        """
        s = r - sigma
        g = epsilon * (r - sigma)
        return g
	
    def eandg(self, i_spin = 1):
        """
        Write a function that computes the Lennard-Jones potential
        V = 0.5 * epsilon * (r-sigma)^2
        dist = self.get_dist_mat()
        """
        self.status.zeros()  
        # parameters
        if i_spin == 1:
            par = self.params[0]
        elif i_spin == 3:
            par = self.params[1]
        else:
            print "no such spin condition!!!"
            exit(1)
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        # print dist
        # sigma epsilon pair
        sigma = par['sigma']; epsilon = par['epsilon'];
        for i in xrange(n_site):
            # force act on site i. f_r dr_i
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos -ibody.pos
                r = np.linalg.norm(rij)
                drdi = - rij / r
                drdj = - drdi
                self.status.potential_energy += self.energy(r, sigma, epsilon)    
                ibody.force += - self.gradient(r, sigma, epsilon) * drdi
                jbody.force += - self.gradient(r, sigma, epsilon) * drdj
        for ibody in sites:
            print "na force", ibody.force
        return
        
    def gdiff(self):
        """ gradient difference of two state """
        self.status.zeros()  
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        # print dist
        # sigma epsilon pair
        #
        g = [np.zeros(3) for i in xrange(n_site)]
        for i in xrange(n_site):
            # force act on site i. f_r dr_i
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos -ibody.pos
                r = np.linalg.norm(rij)
                drdi = - rij / r
                drdj = - drdi
                # spin 1
                par = self.params[0]
                sigma = par['sigma']; epsilon = par['epsilon'];
                t = self.gradient(r, sigma, epsilon) * drdi
                g[i] += t; g[j] += -t;
                # spin 3
                par = self.params[1]
                sigma = par['sigma']; epsilon = par['epsilon'];
                t = - self.gradient(r, sigma, epsilon) * drdi
                g[i] -= t; g[j] -= -t;
        self.status.have(keyword="gdiff", value=g)        
        for v in g:
            print "gdiff", v
        return g
        
    def prob_isc(self, Hso):
        """
        Landau Zener model, ST model
        """
        # gradient
        gdiff = self.status.give(keyword="gdiff")
        # velocity
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        vel = [np.zeros(3) for i in xrange(n_site)]
        for i in xrange(n_site):
            vel[i] = sites[i].vel
        # calc. prob. 
        t = 0.0
        for i in xrange(n_site):
            t += np.dot(gdiff[i], vel[i])
            
        if abs(t) < 1.0e-10:
            prob = 0.0
        else:    
            xi = 8.0 * Hso * Hso / t
            prob = 1.0 - np.exp(-np.pi / 4.0 * xi)
        print "probility, ", prob    
        return prob
        
    def hopping(self, prob):
        ran = np.random.rand()
        i_spin = 1
        if prob > ran:
            i_spin = 3
        print "i_spin, ", i_spin, ran, prob  
        return i_spin
        
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
    status.velocitize(temperature=300)
    ff = ffTwoState(status)
    ff.eandg()
    g = ff.gdiff()
    prob = ff.prob_isc(0.4)
    ff.hopping(prob)    
        
        
        