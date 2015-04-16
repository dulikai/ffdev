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

    def singlePage(filename="abc.xyz"):
        """
        for a xyz file with only one xyz structure
        """
        fp = open(filename, "r")
        n_site = int(fp.readline())
        title = fp.readline().strip()
        sites = [[]for i in xrange(n_site)]
        for i in xrange(n_site):
            line = fp.readline()
            name, pos = self.line2vars(line)
            sites[i].name = name
            sites[i].pos = pos
        fp.close()
        # store vars
        self.sites = sites
        self.n_site = n_site
        fp.close()

          

class bInterpolation():
    def __init__(self):
        
        return

    
    def read(filename=[]):

    
