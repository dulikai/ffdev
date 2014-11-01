#! /usr/bin/env python

#
import os
import numpy as np
from xyzModel import xyzSingle, xyzModel
from bRotation import bRotation as rot


# this class aims to deal with one molecular structure..
class bMole:
    """ molecular class """    
    
    def __init__(self):
        # self.mol = xyz
        
        return
    
    def extend(self):
        pass
        
    def read(self, filename="interface.xyz"):
        """ read how many mol in the xyz file"""
        self.mol = xyzSingle()
        self.mol.read(filename)
        return

    def make_plane(self, mylist):
        geom = self.mol.give(keyword="sites")
        i = int(mylist[0])
        j = int(mylist[1])
        k = int(mylist[2])
        plane = rot.make_plane(geom[i].pos, geom[j].pos, geom[k].pos)
        return plane
        
        
    def rotate(self, pref, pmol):
        """
        get reference plane pref
        get mol plane pmol
        """
        axis, theta = rot.get_axis_theta(pref, pmol)
        print axis, theta
        sites = self.mol.give(keyword="sites")
        # rot
        for ibody in sites:
            pos = ibody.pos
            ibody.pos = rot.rotate_ar_vector(pos, axis, theta)
        return

    def rotate(self, pref, pmol):
        """
        get reference plane pref
        get mol plane pmol
        """
        axis, theta = rot.get_axis_theta(pref, pmol)
        m = rot.get_mat4v(pref, pmol, np.array([0,3,0]))
        print axis, theta
        sites = self.mol.give(keyword="sites")
        # rot
        # for ibody in sites:
            # pos = ibody.pos
            # pos = np.append(pos, 1.0)
            # ibody.pos = rot.rotate_ar_vector(pos, axis, theta) 
            # ibody.pos = np.dot(m, pos)[0:3]
        self.mol.center_id = -1
        self.mol.transform(m)    
        return
        
        
    def dump(self, filename="mol.xyz"):
        """ file name """
        sites = self.mol.give(keyword="sites")
        n_site = self.mol.give(keyword="n_site")
        fp = open(filename, "w")
        print >>fp, "%-10d" % n_site
        print >>fp, ""
        for ibody in sites:
            pos = ibody.pos
            name = ibody.name
            print >>fp, "%-5s%12.6f%12.6f%12.6f" % (name, pos[0], pos[1], pos[2])
        fp.close()    
        return
        
if __name__ == "__main__":

    ref_a = np.array([0.0, 0.0, 0.0])
    ref_b = np.array([1.0, 0.0, 0.0])
    ref_c = np.array([0.0, 0.0, 1.0])
    pref = rot.make_plane(ref_a, ref_b, ref_c)  
    
    # pref = [3, 0, 0, 0] 
    line = input("enter working directory: \n > ")
    print line
    os.chdir(line)
    
    line = input("enter the xyz file name: \n > ")
    xyzfile = line.strip()
    line = input("enter three atom to define molecule plane: \n > ")
    pa, pb, pc = line
    #
    b = bMole()
    b.read(filename=xyzfile)
    pmol = b.make_plane([pa, pb, pc])
    b.rotate(pref, pmol)
    b.dump()


        
        