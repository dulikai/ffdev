#! /usr/bin/env python

#
import os
import numpy as np
from xyzModel import xyzSingle
from bRotation import bRotation as rot


# this class aims to deal with one molecular structure..
class bMole:
    """ molecular class """    
    
    def __init__(self):
         
        return
    
    def extend(self):
        pass
    def read(self, filename="phbr.xyz"):
        """ read how many mol in the xyz file"""
        self.mol = xyzSingle()
        self.mol.read(filename)
        return

    def make_plane(self, i, j, k):
        geom = self.mol.sites
        pi = geom[i].pos
        pj = geom[j].pos
        pk = geom[k].pos
        plane = rot.make_plane(pi, pj, pk)
        return plane
        
        
    def rotate(self, pref, pmol):
        """
        get reference plane pref
        get mol plane pmol
        """
        axis, theta = rot.get_axis_theta(pref, pmol)
        sites = self.mol.give(keyword="sites")
        # rot
        for ibody in sites:
            pos = ibody.pos
            ibody.pos = rot.rotate_ar_vector(pos, axis, theta)
        return
        
    def dump(self, filename="mol.xyz"):
        """ file name """
        sites = self.mol.give(keyword="sites")
        n_site = self.mol.give(keyword="n_site")
        fp = open("mol.xyz", "w")
        print >>fp, "%-10d" % n_site
        print >>fp, ""
        for ibody in sites:
            pos = ibody.pos
            name = ibody.name
            print >>fp, "%-5s%12.6f%12.6f%12.6f" % (name, pos[0], pos[1], pos[2])
        fp.close()    
        return
        
if __name__ == "__main__":
    ilist = [1]
    jlist = [2]
    klist = [3]
    mylist = [ilist, jlist, klist]

    pi = np.array([0.0, 0.0, 0.0])
    pj = np.array([1.0, 0.0, 0.0])
    pk = np.array([0.0, 0.0, 1.0])
    pref = rot.make_plane(pi, pj, pk)
    
    #pref = [3, 0, 0, 0]    

    xyzfile = "phbr.xyz"

    line = input("enter working directory: \n > ")
    os.chdir(line)

    pa, pb, pc = input("enter list to define molecular plane: \n > ")
    b = bMole()
    b.read(filename=xyzfile)
    
    pmol = b.make_plane(pa, pb, pc)
    
    b.rotate(pref, pmol)
    b.dump(filename="mol.xyz")


        
        
