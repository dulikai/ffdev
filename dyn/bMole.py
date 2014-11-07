#! /usr/bin/env python

#
import os
import numpy as np
from xyzModel import xyzSingle
from bRotation import bRotation as rot


class atomicPoint:
    """Hold a 3 dimensional Cartesian coordinate vector."""
    name = "Ar"
    pos = np.zeros(3)

    def __init__(self, coord=np.zeros(3), name="C"):
        """Initialise coordinates."""
        self.pos = coord[0:3]
        self.name = name
        self.frag_id = 0

    def __repr__(self):
        """Print complete string representation of vector."""
        return 'XYZPoint(' + repr(self.pos[0]) + ', ' + repr(self.pos[1]) + ', ' + \
               repr(self.pos[2]) + ')'

    def __str__(self):
        """Return simple string representation of vector to 5dp."""
        return '(%.6f, %.6f, %.6f)' % (self.pos[0], self.pos[1], self.pos[2])

    def __add__(self, other):
        """Overloads the '+' operator to add two vectors."""
        pos = self.pos + other.pos
        return xyzPoint(pos, self.name)

    def __iadd__(self, other):
        """Overloads the '+=' operator to add other to self in place."""
        self.pos += other.pos
        return
        
    def __sub__(self, other):
        """Overloads the '-' operator to subtract two vectors."""
        pos = self.pos - other.pos
        return xyzPoint(pos, self.name)

    def __isub__(self, other):
        """Overloads the '-=' operator to subtract other from self in place."""
        self.pos -= other.pos
        return

    def scale(self, factor):
        """Return the vector multiplied by a scalar factor."""
        factor = float(factor)
        pos = self.pos * factor
        return xyzPoint(pos, self.name)

    def dotProduct(self, other):
        """Returns the dot (scalar) product of two vectors."""
        dotp = np.dot(self.pos, other.pos)
        return xyzPoint(dotp, self.name)

    def crossProduct(self, other):
        """Returns the cross (vector) product of two vectors."""
        pos = np.cross(self.pos, other.pos)
        return xyzPoint(pos, self.name)   

    def getNorm(self):
        """Return the euclidean norm of the vector."""
        norm = np.linalg.norm(self.pos)
        return norm

    def normalise(self):
        """Return the vector after normalisation.

        If the normalisation factor is zero, the original (zero-)vector is
        returned unchanged.
        """
        try:
            normFactor = 1.0 / self.getNorm()
        except ZeroDivisionError:
            normFactor = 1.0
        return self.scale(normFactor)
        
    def transfrom(self, m):
        """ return m x point
        transformation of a point
        """
        pos = np.append(self.pos, 1.0)
        pos = np.dot(m, pos)
        return xyzPoint(pos[0:3], self.name)
            
            

# this class aims to deal with one molecular structure..
class bMole:
    """ molecular class """    
    
    def __init__(self):
         
        return
    
    def slice(self, sitelist):
        """ a list of mole """
        newmol = bMole()
        for i in sitelist:
            mol.append()
    
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


        
        
