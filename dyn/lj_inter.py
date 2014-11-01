#! /usr/bin/env python

import numpy as np
import random
# The following parameters will be used for the Argon atom:
 # mass : 39.948 g/mol
# epsilon = 0.0661 j/mol
# sigma = 0.3345 nm
# The atom coordinates will be expressed in nm.


class lj_inter:
    def __init__(self):
        self.cluster = {}
        self.dist_mat = []
        
        return

        
    def potential(self, sigma, epsilon):
        """
        Write a function that computes the Lennard-Jones potential
        V = 4 \epsilon [(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6]
        """
        lj.get_distances()
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
        
    def kinetic(self):
        """
        T = 0.5 \sum {m v^2}
        """
        n_site = self.cluster['n_site']
        site = self.cluster['site']
        t = 0.0
        for i in xrange(n_site):
            vel = site[i]['vel']
            mass = site[i]['mass']
            t += mass * np.dot(vel, vel) * 0.5
        return t    
       
    def interaction(self, sigma, epsilon):
        """
        whole system
        """
        vp = self.potential(sigma, epsilon)
        t = self.kinetic()
        v = vp + t
        return v
        
    def modify_site_coord(self, isite, idim, delta):
        """
        for numerical diff; modify a var
        """
        self.cluster['site'][isite]['coord'][idim] += delta
        return
        
    def compute_nd_acc(self, sigma, epsilon):
        """
        numerical force
        """
        delta = 1.0e-6
        n_site = self.cluster['n_site']
        site = self.cluster['site']
        for i_site in xrange(n_site):
            mass = site[i_site]['mass']
            acc = np.array([0.0, 0.0, 0.0])
            for j_dim in xrange(3):
                self.modify_site_coord(i_site, j_dim, delta)
                v2 = self.potential(sigma, epsilon)
                self.modify_site_coord(i_site, j_dim, -2*delta)
                v1 = self.potential(sigma, epsilon)
                self.modify_site_coord(i_site, j_dim, delta)
                acc[j_dim] = - (v2 - v1) / (2.0 * delta) / mass
            site[i_site]['acc'] = acc
            # print "acc",acc
        return
            
    def read_site_coord(self):
        """
        read in coordinate of site
        """
        fp = open("site.xyz", "r")
        line = fp.readline()
        n_atom = int(line.strip())
        line = fp.readline()
        site = []
        for i in xrange(n_atom):
            line = fp.readline()
            # print line
            rec = line.split()
            name = rec[0]
            coord = np.array([float(rec[1]), float(rec[2]), float(rec[3])])
            vel = np.array([0.0, 0.0, 0.0])
            site.append({'name': name, 'coord': coord, 'vel': vel})
        self.cluster = {'site': site, 'n_site': n_atom}
        return
        
    def read_site_mass(self):
        """
        read in site mass from a file
        """
        cluster = self.cluster
        fp = open("site.mass", "r")
        line = fp.readline()
        n_site = int(line.strip())
        if cluster['n_site'] != n_site:
            print "error: check n_site consitancy"
            exit(1)
        line = fp.readline()
        for i in xrange(n_site):
            isite = cluster['site'][i]
            line = fp.readline()
            mass = float(line)
            isite['mass'] = mass
        return
            
    def get_distances(self):
        """
        calculate the dist matrix.
        """
        site = self.cluster['site']
        n_site = self.cluster['n_site']
        self.dist_mat = [[0.0 for j in xrange(n_site)] for i in xrange(n_site)]
        for i in xrange(n_site):
            for j in xrange(n_site):
                coordi = site[i]['coord']
                coordj = site[j]['coord']
                self.dist_mat[i][j] = np.linalg.norm(coordi-coordj)
        
        return
        
    def verlet(self, sigma, epsilon):
        dt = 0.1
        n_site = self.cluster['n_site']
        site = self.cluster['site']
        energy = self.interaction(sigma, epsilon)
        self.compute_nd_acc(sigma, epsilon)
        fp = open("output.xyz", "w")
        for i in xrange(1000):
            energy = self.interaction(sigma, epsilon)
            self.compute_nd_acc(sigma, epsilon)
            print >>fp, "%10d" % n_site
            print >>fp, "Step %d E=%20.12f" % (i, energy)
            print "Step %d E=%20.12f" % (i, energy)
            acc_old = []
            for i_site in xrange(n_site):
                vel = site[i_site]['vel']
                acc = site[i_site]['acc']
                coord = site[i_site]['coord']
                print >>fp, "Ar %12.6lf%12.6lf%12.6lf" % (coord[0]*10, coord[1]*10, coord[2]*10) 
                coord += vel * dt + acc * dt * dt * 0.5
                # print site[i_site]['coord']            
                acc_old.append(acc.copy())
                
            self.compute_nd_acc(sigma, epsilon)
            for i_site in xrange(n_site):
                vel = site[i_site]['vel']
                acc = site[i_site]['acc']
                vel += (acc_old[i_site] + acc) * dt * 0.5 
                # print acc_old-acc

        return
if __name__ == "__main__":
    lj = lj_inter()
    sigma = 0.3345
    epsilon = 0.0661
    lj.read_site_coord()
    lj.read_site_mass()
    # line = raw_input("enter inter-atom distances ? \n > ")
    # r = float(line)
    # v = lj.interaction(sigma, epsilon)
    # lj.compute_nd_acc(sigma, epsilon)
    # print lj.dist_mat
    lj.verlet(sigma, epsilon)
    
    # print v
    