#! /usr/bin/env python

import os
import copy
import numpy as np
from bStatus import bStatus
from bRotation import bRotation
from xyzModel import xyzSingle, xyzModel

# symmetry builder

class bPolygon:
    status = bStatus()
    polygon = {}
    
    def regular_polygon(self, n=3, rad=3.0, theta=0.0, theta0=0.0, mode="polygon"):
        """ build regular polygonal distribution 
            job may 'circle' or 'polygen'
        """
        if n < 3:
            print "n < 3 in regular_polygon (??) !!!"
            exit(1)
        theta = theta / 180.0 * np.pi
        theta0 = theta0 / 180.0 * np.pi
        dt = 2.0 * np.pi / n
        # dist = np.sqrt(2.0 * rad * rad * (1 - np.cos(dt)))
        # print dist
        if mode == "polygon":
            rad = rad / np.sqrt(2.0*(1-np.cos(dt)))
            print rad
        # polygon vertices
        vertices = []
        for i in xrange(n):
            t = theta0 + i * dt
            x = rad * np.cos(t)
            y = rad * np.sin(t)
            z = 0.0
            coord = np.array([x, y, z])
            vertices.append(coord)
        # polygon edges
        edges = []
        for i in xrange(n):
            j = i + 1
            if i == n - 1:
                j = 0
            edges.append(vertices[i] - vertices[j])
            # print vertices[i], vertices[j], i,j
        # the branches, that is X ---> A vector
        branches = []
        axis = np.array([0., 0., 1.])   # norm of the polygon plane (z axis)
        for i in xrange(n):
            vec = edges[i] 
            v = bRotation.rotate_ar_vector(vec, axis, theta) 
            print theta
            # print v, vec, axis, theta
            branches.append(v)
        plane = bRotation.make_plane(vertices[0], vertices[1], vertices[2])
        self.polygon = {'n_sides': n, 'vertices': vertices, 'edges': edges, 
                        'branches': branches, 'plane': plane}
        return
        
    def mapping(self, origin, vec):
        vertices = self.polygon['vertices']
        branches = self.polygon['branches']
        n_sides = self.polygon['n_sides']
        mat = []
        for i in xrange(n_sides):
            # print vec, branches[i], origin, vertices[i]
            a, b = bRotation.get_axis_theta(vec, branches[i])
            va = np.array([1.0, 1.0, 0.0])
            vb = np.array([1.0, 0.0, 0.0])
            m = bRotation.get_mat4t(vec, branches[i], origin, vertices[i]) 
            m = bRotation.get_mat4t(vec, branches[i], origin, origin) 
            # m = bRotation.get_mat4t(va, vb, origin, vertices[i]) 
           
            mat.append(m)
        self.polygon['mapping'] = mat
        return
        
    def dump(self, filename="polygon.xyz"):
        fp = open(filename, "w")
        vertices = self.polygon['vertices']
        branches = self.polygon['branches']
        plane = self.polygon['plane']
        n = self.polygon['n_sides']
        print >>fp, "%-10d" % (n*2)
        print >>fp, ""
        for i in xrange(n):
            v = vertices[i]
            b = branches[i] + v
            print >>fp, "%-5s%12.6f%12.6f%12.6f" % ("I", v[0], v[1], v[2])
            print >>fp, "%-5s%12.6f%12.6f%12.6f" % ("C", b[0], b[1], b[2])
        print >>fp, "PLANE"    
        for i in xrange(4):
            print >>fp, "%12.6f" % (plane[i]),
        fp.close()
        return
    
    def images(self, xyz):
        mat = self.polygon['mapping']
        mol = xyzSingle()
        for m in mat:
            t = copy.deepcopy(xyz)
            t.transform(m)
            # t.dump(filename="tmp.xyz")
            mol.extend(t)
        return mol
    
    
    def build_test(self, config):
        # mole
        core = config['mole']['core']
        direction = config['mole']['direction']
        
        # shape
        n_sides = config['polygen']['n_sides']
        t0 = config['polygen']['theta0']
        mode = config['polygen']['mode']
        # min max step
        rmin = config['radii']['min']
        rmax = config['radii']['max']
        rsize = config['radii']['size']
        tmin = config['theta']['min']
        tmax = config['theta']['max']
        tsize = config['theta']['size']
        # files
        xyzfile = config['file']['xyz']
        #
        xyz = xyzSingle()
        xyz.read(filename=xyzfile)
        origin, vec = xyz.set_info(origin=core, direction=direction)
        #
        rnum = int((rmax-rmin) / rsize)
        tnum = int((tmax-tmin) / tsize)  
        # 
        polygon = bPolygon()
        model = xyzModel()
        v = np.array([0.0, 1.0, 0.0])
        for j in xrange(tnum):
            t = tmin + j * tsize 
            mx = bRotation.get_rot_mat4v(v, vec)
            polygon.mapping(origin, vec)

                # polygon.print_polygon(str(i)+str(j)+".xyz")
            newxyz = polygon.images(xyz)
            model.extend(newxyz)
                # print len(newxyz.sites)
        model.dump()
        # print model.n_xyz
        return    
    def build(self, config):
    
        # mole
        core = config['mole']['core']
        direction = config['mole']['direction']
        
        # shape
        n_sides = config['polygen']['n_sides']
        t0 = config['polygen']['theta0']
        mode = config['polygen']['mode']
        # min max step
        rmin = config['radii']['min']
        rmax = config['radii']['max']
        rsize = config['radii']['size']
        tmin = config['theta']['min']
        tmax = config['theta']['max']
        tsize = config['theta']['size']
        # files
        xyzfile = config['file']['xyz']
        #
        xyz = xyzSingle()
        xyz.read(filename=xyzfile)
        origin, vec = xyz.set_info(origin=core, direction=direction)
        #
        rnum = int((rmax-rmin) / rsize)
        tnum = int((tmax-tmin) / tsize)  
        # 
        polygon = bPolygon()
        model = xyzModel()
        for i in xrange(rnum):
            for j in xrange(tnum):
                r = rmin + i * rsize
                t = tmin + j * tsize        
                polygon.regular_polygon(n=n_sides, rad=r, theta=t, theta0=t0, mode=mode)
                polygon.mapping(origin, vec)

                # polygon.print_polygon(str(i)+str(j)+".xyz")
                newxyz = polygon.images(xyz)
                model.extend(newxyz)
                # print len(newxyz.sites)
        model.dump()
        # print model.n_xyz
        return
        
#        
# the rotation is anti-clockwise for positive angle 
# the axis may be zero vector if the vec and branches 
# are in the same line
# which is possible in some case, if this is o, change polygon-theta0..
#
#  positive angle: rotate from right to left
#  define the default direction as:
#
#          /   
#         /
#        /\
#       /  \
#      /    \
#  ----------\
#             \
#
#        
if __name__ == "__main__":

    line = input("enter working directory: \n > ")
    print line
    os.chdir(line.strip())
    
    # template
    # xyz = xyzSingle()
    # xyz.read(filename="mol.xyz")
    # xyz.dump(filename="tmp.xyz")
    # origin, vec = xyz.set_info(origin=11, direction=5)
    
    # polygon = bPolygon()
    # polygon.regular_polygon(n=3, rad=8.0, theta=30.0, theta0=0.0)
    # polygon.dump(filename="polygon.xyz")
    
    # polygon.mapping(origin, vec)
    
    # newxyz = polygon.images(xyz)
    
    # newxyz.dump(filename="newxyz.xyz")
    
    
    # exit(1)
   
    # polygon.print_polygon()
   
    
    
    config = {}
    config['radii'] = {'min': 3, 'max': 6, 'size': 0.1}
    config['theta'] = {'min': 0, 'max': 120, 'size': 10}
    config['polygen'] = {'n_sides': 3, 'theta0': 8, 'mode': 'circle'}
    config['mole'] = {'core': 11, 'direction': 5}
    config['file'] = {'xyz': "mol.xyz"}
    polygon = bPolygon()
    polygon.build(config)
    
    
    
    
    # origin = np.array([0, 1, 1])
    # vec = np.array([3,3,3])
    # b = bPolygon()
    # b.regular_polygon(n=3, rad=3.0, theta=0.0, theta0=0.0)
    # b.mapping(origin, vec)

    
    