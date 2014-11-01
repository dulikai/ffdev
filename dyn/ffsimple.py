#! /usr/bin/env python

# force field model: simple case

import numpy as np

from bSite import bSite

class ffsimple():
    """
    a very simple case of particle interaction
    in which,
    the trans. freedom is interaction with l
    """
    def __init__(self, bSite):
        
        return
        
    def lj_pot(self, epsilon, sigma):
        """
        V_{rep-disp} = 4 \epsilon \cdot [V_{rep} - V_{disp}]
        """
        
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
            print line
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
        
    # trans. and rot.    
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