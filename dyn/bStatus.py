#! /usr/bin/env python

import os
import numpy as np

#
# site layer page node cluster group 
#

class bStatus:
    """
    layer is composed by sites
    """
    box = []
    sites = []
    n_site = 0
    n_freedom = 0
    potential_energy = 0.0
    
    def __dir__(self):
        return ["sites", "n_site"]
        
    def __init__(self, config = {}):
        """ read in restart file """
        self.read_box()
        self.read_coord()
        self.read_mass()
        return
        
    def give(self, keyword=""):
        try:
            value = getattr(self, keyword)
        except AttributeError:
            print("@error: no such model keyword (bStatus)!!!")
            exit(1)
        return value
        
    def have(self, keyword="", value=""):
        setattr(self, keyword, value)
        return

    def zeros(self):
        self.potential_energy = 0.0
        sites = self.give(keyword="sites")
        for ibody in sites:
            ibody.force = np.zeros(3)
            ibody.torque = np.zeros(3)
        return
        
    def dump(self, fp):
        sites = self.sites  
        n_site = self.n_site
        print >>fp, "%10d" % n_site
        print >>fp, ""
        # dump pos & vel
        for ibody in sites:
            name = ibody.name
            pos = ibody.pos
            ang = ibody.angpos
            vel = ibody.vel
            acc = ibody.force
            print >>fp, "%10s%12.6f%12.6f%12.6f" % (name, pos[0], pos[1], pos[2]),
            print >>fp, "%12.6f%12.6f%12.6f" % (ang[0], ang[1], ang[2])
            # print >>fp, "%10s%12.6f%12.6f%12.6f" % (name, acc[0], acc[1], acc[2])
    
        return
        
    @staticmethod    
    def line2vars(line):
        rec = line.split()
        vel = np.array([0.0, 0.0, 0.0])
        name = rec[0]
        coord = np.array([float(rec[1]), float(rec[2]), float(rec[3])])
        if len(rec) == 7:
            vel = np.array([float(rec[4]), float(rec[5]), float(rec[6])])
        return name, coord, vel
        
    def read_box(self):
        """
        read in pbc condition
        """
        box = []
        try:
            os.path.isfile("box.pbc")
            fp = open("box.pbc", "r")
            line = fp.readline()
            rec = line.split()
            for val in rec:
                box.append(float(val))
        except IOError:
            print "NO PBC FILE"
        self.box = np.array(box)    
        return
                
    def read_coord(self):
        """
        read cart. coordinate in xyz format
        read euler angle in abc format
        """
        # efp xyz
        fp = open("efp.xyz", "r")  
        n_site = int(fp.readline())
        line = fp.readline()
        sites = [bSite() for i in xrange(n_site)]
        for i in xrange(n_site):
            line = fp.readline()
            name, pos, vel = self.line2vars(line)  
            sites[i].name = name
            sites[i].pos = pos
            sites[i].vel = vel   
        fp.close()
        # efp abc
        fp = open("efp.abc", "r")
        n_site = int(fp.readline())
        line = fp.readline()
        for i in xrange(n_site):
            line = fp.readline()
            name, angpos, angvel = self.line2vars(line) 
            sites[i].angpos = angpos
            sites[i].angvel = angvel
        # store vars
        self.sites = sites
        self.n_site = n_site
        self.n_freedom = n_site * 6
        fp.close()
        # done
        return
        
    def read_mass(self):
        """
        read in site mass
        """
        sites = self.sites
        # read mass
        fp = open("efp.mass", "r")
        n_site = int(fp.readline())
        line = fp.readline()
        for i in xrange(n_site):
            mass = float(fp.readline())
            sites[i].mass = mass 
        fp.close()        
        # read inertia
        fp = open("efp.inertia", "r")
        n_site = int(fp.readline())
        line = fp.readline()
        for i in xrange(n_site):
            rec = fp.readline().split()
            inertia = np.array([float(rec[0]), float(rec[1]), float(rec[2])])
            sites[i].inertia = inertia
        fp.close()    
        return
    
    def print_info(self):
        sites = self.sites
        for mysite in sites:
            print mysite.pos, mysite.angpos
        
class bSite:
    """
    typical structure of rigid body like site model
    every site can be a rigid body
    """
    name = "C"
    # center of mass: trans ...
    mass = 0.0
    pos = np.zeros(3)
    vel = np.zeros(3)
    force = np.zeros(3)
    acc = np.zeros(3)    # acceleration
    #rotation ...
    inertia = np.zeros(3)
    angpos = np.zeros(3)
    angvel = np.zeros(3)
    torque = np.zeros(3)
    angacc = np.zeros(3)
    # angpos_old = np.zeros(3)
    # angvel_old = np.zeros(3)
    inertia_inv = np.zeros(3)    
    # angular momentum
    rot = np.zeros(3)
    angmom = np.zeros(3)
    # angmom_old = np.zeros(3)
    
    def __str__(self):
        mystr = "%32s" % self.name
        return mystr 
    
    def dump_xyz(self):
        pos = self.pos
        mystr = "%20.8f%20.8f%20.8f\n" % (pos[0], pos[1], pos[2])
        vel = self.vel
        mystr += "%20.8f%20.8f%20.8f\n" % (vel[0], vel[1], vel[2])
        angpos = self.angpos
        mystr += "%20.8f%20.8f%20.8f\n" % (angpos[0], angpos[1], angpos[2])
        angvel = self.angvel
        mystr += "%20.8f%20.8f%20.8f\n" % (angvel[0], angvel[1], angvel[2])
        return
        
    def have(self, keyword="", value=""):
        setattr(self, keyword, value)
        return
    
if __name__ == "__main__":
    site = bSite()
    print site.mass
    status = bStatus()
    print 
    sites = status.give(keyword = "sites")
    print sites[2]
    # status.print_info()

    
    
    
    
    

        