#! /usr/bin/env python

import numpy as np
from xyzSingle import xyzSingle
from bPolygon import bPolygon



class bBuilder:
    """ build the images """
    
    
    def __init__(self):
        xyz = xyzSingle()
        xyz.read()
        origin = xyz.point(1)
        vec = xyz.bond_vector(1,2)        
        polygon = bPolygon()
        polygon.regular_polygon(n=3, rad=3.0, theta=0.0, theta0=0.0)
        polygon.mapping(origin, vec)
        self.images(polygon, xyz)
        
        


        
    def images(self, polygon, xyz):
        mats = polygon.polygon['mapping']
        mol = xyzSingle()
        for m in mats:
            template = copy.deepcopy(xyz)
            t = template.transfrom(m)
            mol.extend(mol)
        return mol
        
if __name__ == "__main__":
    b = bBuilder()
    
