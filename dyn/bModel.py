#! /usr/bin/env python

# this class use point representation of reactive site..
#
import numpy as np

# point for a single reactive site.
class bPoint():
    label = "UNK"
    coord = np.zeros(3)
    mass = -1.0
    charge = 0.0
    epsilon = 0.0
    sigma = 0.0
    
class bFragment():
    # fragment name
    name = "UNK"
    # fragment center of mass
    coord = np.zeros(3)
    # rotation matrix representing orientation of a fragment
    rotmat = np.zeros((3,3))
    # rotation euler angle representing orientation of a fragment
    euler = np.zeros(3)
    # number of points in this fragment
    n_points = 0
    # fragment points
    points = []
    # 
    def setup(self):
        self.points = [bPoint() for i in xrange(self.n_points)]
        
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
                
    
    