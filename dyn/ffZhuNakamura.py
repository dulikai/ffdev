#! /usr/bin/env python
#
# ----------------------------------------------------
# @ dulikai
# @ 2014.11.25
# @ Zhu-Nakamura formula
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
# suppose direction is observed.

# find beta

# estimate minimum

# calc. alpha

# calc. PZN


class ffTwoState:
    """
    a very simple case of particle interaction
    in which,
    harmonic potential
    """
    def __init__(self, status):
        self.status = status
        # sigma epsilon pairs
        self.params = {}
        self.params['S_0'] = {'sigma': 2.345, 'epsilon': 1.661}   # spin singlet
        self.params['T_1'] = {'sigma': 3.500, 'epsilon': 2.500}  # spin triplet
        return

    @staticmethod
    def func(x):
        return x*x - 3.0
    @staticmethod    
    def dfunc(x):
        df = 2*x
        if abs(df) <= 1.0e-10:
            df = 1.0
        return df
        
    def find_root(self, x0=100):
        x1 = x0 - self.func(x0)/self.dfunc(x0)
        eps = 1.0e-8
        print x0, x1
        while abs(x1-x0) >= eps or self.func(x0) >= eps:
            # raw_input("contineu..")
            print "x0: ", x0
            x0 = x1
            x1 = x0 -  self.func(x0)/self.dfunc(x0) 
            print "x1: ", x1
            # print x0, x1, abs(x0-x1), self.func(x0), self.dfunc(x0), self.func(x0)/self.dfunc(x0)
        return  x1
        
        
    def detect_min(self, ene1, ene2):
        """
        vh: previous pes; v: current pes; vl: next pes point value
        """
        de = ene1 - ene2
        v = np.abs(de)
        flag = "NO"
        if v[1] < v[0] and v[1] < v[2]:
            flag = "YES"
        return flag    
        
        
        
    def get_direction(self):
        grad_diff(self, i_spin="S_0", j_spin="T_1"):  

        
    def energy_diff(self, i_spin="S_0", j_spin="T_1"):
        """
        energy difference
        """        
        i_spin = self.get_curr_state()
        # parameters
        if i_spin in self.params.keys():
            par = self.params[i_spin]
        else:
            print "no such spin condition!!!"
            exit(1)
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        self.status.potential_energy = 0.0
        # print dist
        # sigma epsilon pair
        ediff = 0.0
        sigma = par['sigma']; epsilon = par['epsilon'];
        for i in xrange(n_site):
            # force act on site i. f_r dr_i
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos -ibody.pos
                r = np.linalg.norm(rij)
                par = self.params[i_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                ediff += self.energy(r, sigma, epsilon) 
                par = self.params[j_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                ediff -= self.energy(r, sigma, epsilon) 
        return ediff
            
    
    def energy_twostate(self, alpha, i_spin="S_0", j_spin="T_1"):
        """
        energy difference
        """        
        i_spin = self.get_curr_state()
        # parameters
        if i_spin in self.params.keys():
            par = self.params[i_spin]
        else:
            print "no such spin condition!!!"
            exit(1)
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        direction = self.status.give(keyword="direction") 
        direction *= alpha
        self.status.potential_energy = 0.0
        # print dist
        # sigma epsilon pair
        ediff = 0.0
        for i in xrange(n_site):
            # force act on site i. f_r dr_i
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos - ibody.pos + direction[i] - direction[j]
                r = np.linalg.norm(rij)
                par = self.params[i_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                ediff += self.energy(r, sigma, epsilon) 
                rij = jbody.pos - ibody.pos - direction[i] + direction[j]
                r = np.linalg.norm(rij)
                par = self.params[j_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                ediff -= self.energy(r, sigma, epsilon) 
        return ediff
            
    def grad_twostate(self, alpha, i_spin="S_0", j_spin="T_1"):
        """ 
        gradient difference of two state
        U+(r0+alphaA)-U-(r0-alphaA) = 0
        AU+', -AU-'
        """
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        direction = self.status.give(keyword="direction") 
        direction *= alpha
        # print dist
        galpha = 0.0
        # sigma epsilon pair
        ga = [np.zeros(3) for i in xrange(n_site)]
        gb = [np.zeros(3) for i in xrange(n_site)]
        for i in xrange(n_site):
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos -ibody.pos + direction[i] - direction[j]
                r = np.linalg.norm(rij)
                drdi = - rij / r
                drdj = - drdi
                par = self.params[i_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                t = self.gradient(r, sigma, epsilon) * drdi
                ga[i] += t; ga[j] += -t;
                rij = jbody.pos - ibody.pos - direction[i] + direction[j]
                r = np.linalg.norm(rij)
                drdi = - rij / r
                drdj = - drdi
                par = self.params[j_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                t = - self.gradient(r, sigma, epsilon) * drdi
                gb[i] -= t; gb[j] -= -t;
        
        for i in xrange(n_state):
            galpha += np.dot(ga[i], direction[i])
            galpha += np.dot(gb[i], direction[i])

        # for v in g:
            # print "gdiff", v
        return galpha
                    
    def get_alpha(self):
        gdiff = self.give(keyword="gdiff")
        vec = gdiff / np.linalg.norm(gdiff)
        alhpa = 1.0
        # f(alpha) = fa(r0+alpha A) - fb(r0-alpha A)
        # see JCP, 2001,115, 11036
        find_root()
        return alpha
        
    def get_energy0(self):
        """
        EP = U+ = U-
        """
        ep = # jcp,2001,115, 11036 eq(2.7)
        return ep
        
        
    def get_d2():
        return
    def get_a2():
        return
    def get_b2():
        return
    #    
    # see also Computational and Theoretical Chemistry
    #  Volume 1023, 1 November 2013, Pages 10¨C18
    # Nonadiabatic dynamics study of bridged-azobenzene 
    # by the time-dependent density functional tight-binding method
    #
    
    def prob_zn(self, a, b):
        
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
	
    def eandg(self):
        """
        Write a function that computes the Lennard-Jones potential
        V = 0.5 * epsilon * (r-sigma)^2
        dist = self.get_dist_mat()
        """        
        self.status.zeros()  
        i_spin = self.get_curr_state()
        # parameters
        if i_spin in self.params.keys():
            par = self.params[i_spin]
        else:
            print "no such spin condition!!!"
            exit(1)
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        self.status.potential_energy = 0.0
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
                ibody.force += -self.gradient(r, sigma, epsilon) * drdi
                jbody.force += -self.gradient(r, sigma, epsilon) * drdj
        for ibody in sites:
            print "na force", ibody.force
        return
        
    def grad_diff(self, i_spin="S_0", j_spin="T_1"):
        """ gradient difference of two state """
        # calc.
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        # print dist
        # sigma epsilon pair
        g = [np.zeros(3) for i in xrange(n_site)]
        for i in xrange(n_site):
            for j in xrange(i+1, n_site):
                ibody = sites[i]; jbody = sites[j]
                rij = jbody.pos -ibody.pos
                r = np.linalg.norm(rij)
                drdi = - rij / r
                drdj = - drdi
                par = self.params[i_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                t = self.gradient(r, sigma, epsilon) * drdi
                g[i] += t; g[j] += -t;
                par = self.params[j_spin]
                sigma = par['sigma']; epsilon = par['epsilon'];
                t = - self.gradient(r, sigma, epsilon) * drdi
                g[i] -= t; g[j] -= -t;
        self.status.have(keyword="gdiff", value=g)        
        # for v in g:
            # print "gdiff", v
        return g
        
    def prob_isc(self, Hso):
        """
        Landau Zener model, ST model
        """
        # gradient
        gd = self.status.give(keyword="gdiff")
        # velocity
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        vel = [np.zeros(3) for i in xrange(n_site)]
        for i in xrange(n_site):
            vel[i] = sites[i].vel
        # calc. prob. 
        t = 0.0
        for i in xrange(n_site):
            t += np.dot(gd[i], vel[i])
            
        if abs(t) < 1.0e-10:
            prob = 0.0
        else:    
            xi = 8.0 * Hso * Hso / t
            prob = 1.0 - np.exp(-np.pi / 4.0 * xi)
        print "probility: %12.6f, xi: %12.6f" % (prob, xi)    
        return prob
        
    def get_curr_state(self):
        """ obtain current electronic state """
        i_spin = self.status.give(keyword="i_spin")
        if i_spin is None:
            i_spin = "S_0"
        return i_spin
        
    def get_new_state(self):
        """ obtain current electronic state """
        # ??? change later
        i_spin = self.status.give(keyword="i_spin")
        if i_spin is None:
            i_spin = "S_0"
        elif i_spin == "S_0":
            i_spin = "T_1"
        else:
            i_spin = "S_0"
        return i_spin
        
    def hopping(self):
        i_spin = self.get_curr_state()
        ES = self.compute_energy(i_spin="S_0")
        ET = self.compute_energy(i_spin="T_1")
        de = abs(ES-ET)
        if i_spin == "S_0":
            x_spin = "S_0"; y_spin = "T_1";
        else:
            x_spin = "T_1"; y_spin = "S_0"; 
        self.grad_diff(i_spin=x_spin, j_spin=y_spin)
        ran = 0.0; prob = 0.0
        if de > 0.5:
            print "large delta E is : %s" % de
        else:
            print "small delta E is : %s" % de
            prob = self.prob_isc(0.6)
            ran = np.random.rand()
            if prob > ran:
                i_spin = self.get_new_state()
        self.status.have(keyword="i_spin", value=i_spin)        
        print "i_spin: %s, ran: %12.6f, prob: %12.6f" % (i_spin, ran, prob)

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
        
    def compute_energy(self, i_spin="S_0"):
        """ COMPUTE INTERACTION ENERGY """
        if i_spin in self.params.keys():
            par = self.params[i_spin]
        else:
            print "no such spin condition!!!"
            exit(1)
        sigma = par['sigma']; epsilon = par['epsilon'];
        #
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
                ene += self.energy(r, sigma, epsilon)    
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
    # status.velocitize(temperature=300)
    ff = ffTwoState(status)
    # ff.eandg()
    # g = ff.grad_diff()
    # prob = ff.prob_isc(0.4)
    # ff.hopping() 
    
    ff.find_root(x0=0)
    a = np.array([1.0, 0.0, 5.0])
    b = np.array([-2.0, -1.0, -3.0])
    print ff.detect_min(a, b)

    
        