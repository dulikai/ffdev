#! /usr/bin/env python

# linear interpolation within two or several geometries.

import os
import numpy as np


class xyzRoom:
    """
    yet another xyz file parser
    """
    def __init__(self):
        self.vars = {}

        return

    @staticmethod
    def line2vars(line):
        rec = line.split()
        vel = np.array([0.0, 0.0, 0.0])
        name = rec[0]
        coord = np.array([float(rec[1]), float(rec[2]), float(rec[3])])
        if len(rec) == 7:
            vel = np.array([float(rec[4]), float(rec[5]), float(rec[6])])
        return name, coord


          

class bInterpolation():
    def __init__(self):
        
        return

    def line2site(line = "C 0.0 0.0 0.0"):
        record = line.split()
        vel = np.array([0.0, 0.0, 0.0])
        name = record[0]
        coord = [float(f) for f in record[1:4]]
        coord = np.array(coord)

        return name, coord

    def singlexyz(filename = "abc.xyz"):
        """ read in a xyz format file """
        fp = open(filename, "r")
        n_site = int(fp.readline())
        title = fp.readline().strip()
        sites = []
        for i in xrange(n_site):
            line = fp.readline()
            record = line.split()
            name = record[0]
            coord = np.array([float(f) for f in record[1:4]])
            sites.append({'name': name, 'coord': coord})
        fp.close()
             
        return sites

    def get_new_pos(self, pointA, pointB, i, n):
        """
        get the new position via linear interpolation
        ith point from A, 
        """
        k = float(i) / n
        p = pointA * (1 - k) + pointB * k
        return p

    def get_new_geom(self, beg, end, i, n):
        """
        given structural beg & end
        generate n geometries
        """
        n_site = 0
        sites = []
        for xb, xe in zip(beg, end):
            coordb = xb['coord']
            

        return

    
    

    
