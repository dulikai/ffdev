#! /usr/bin/env python

import os
import copy
import numpy as np
from bStatus import bStatus
from bRotation import bRotation
from xyzModel import xyzSingle, xyzModel

# PES SCAN generator for 2d and 1d.
# bond angle dihedral etc
# the molecule is seperated into different fragment.

class bPESGrid:
    # status = bStatus()
    topology = {}
    
    def get_coord(self, idx):
        """
        idx is a list of atoms, return the geom. center
        """
        n = len(idx)
        if n == 0:
            print "empty atom list, error.."
            exit(1)
            
    def radical(self, xyz, tbl={} ):
        """
        group/fragment 2 is the moving group by default
        """
        # suppose the xyz geom is given..
        
        frg = [[0,1,2,3], [4]]
        ndx = [[0], [4]]        
        tbl = {"start": 3.0, "end": 10.0, "size": 0.1}
        
        # build ndx center for radical moving..
        cA = xyz.get_center(ndx[0])
        cB = xyz.get_center(ndx[1])
        
        # mapping   
        # setup directional vector..
        direction = cB - cA
        udir = direction / np.linalg.norm(direction)
        # place frg 1 fixed.
        coordA = cA
        # then we can build a series of A, B pairs,
        # for which, the radical scan is done..
        #
        start = tbl['start']
        end = tbl['end']
        size = tbl['size']
        n_points = int((end-start)/size)
        pairs = []
        for i in xrange(n_points):
            rad = start + size * i
            coordB = cA + udir * rad
            pairs.append([coordA, coordB])
            
        # build transformation matrix
        mat = []
        for i in xrange(n_points):
            coordA = pairs[i][0]
            coordB = pairs[i][1]
            vec = coordA - cA
            mA = bRotation.get_trans_mat4(vec)
            vec =  coordB - cB
            mB = bRotation.get_trans_mat4(vec)
            mat.append([mA, mB])
            
        # do the transformation
        model = xyzModel()
        for m in mat:
            t = copy.deepcopy(xyz)
            t.transform(m[0], frg[0])
            t.transform(m[1], frg[1])
            model.extend(t)
        model.dump()    
        return model
        
        
    def angle(self, mol):
        """
        frag 1 is fixed by default..
        """
         # suppose the xyz geom is given..
        
        frg = [[0,1,2,3,4], [5,6,7,8,9]]
        ndx = [[0], [4], [9]]        
        tbl = {"start": 90.0, "end": 270.0, "size": 1.0}
        
        # build ndx center for angle moving..
        cA = xyz.get_center(ndx[0])
        cB = xyz.get_center(ndx[1])
        cC = xyz.get_center(ndx[2])
        print cA, cB, cC
        # mapping   
        # setup 
        plane = bRotation.make_plane(cA, cB, cC)
        axis = plane[0:3]
        # place frg 1 fixed.
        coordA = cA
        coordB = cB
        # then we can build a series of A, B pairs,
        # for which, the radical scan is done..
        #
        start = tbl['start'] / 180.0 * np.pi
        end = tbl['end'] / 180.0 * np.pi
        size = tbl['size'] / 180.0 * np.pi
        n_points = int((end-start)/size)
        pairs = []
        for i in xrange(n_points):
            theta = start + size * i
            xA = np.array([-np.linalg.norm(cA-cB), 0.0, 0.0])
            xB = np.zeros(3)
            r0 = np.linalg.norm(cB-cC)
            x = r0*np.cos(np.pi-theta)
            y = r0*np.sin(np.pi-theta)
            xC = np.array([x, y, 0.0])
            #            
            xBA = xA - xB
            cBA = cA - cB
            m4 = bRotation.get_mat4t(xBA, cBA, cB, xB)
            coordC = np.dot(m4, np.append(xC, 1.0))[0:3]
            pairs.append([coordA, coordB, coordC])
            print xB, xC
        # exit()    
        # build transformation matrix
        mat = []
        for i in xrange(n_points):
            coordA = pairs[i][0]
            coordB = pairs[i][1]
            coordC = pairs[i][2]
            # A-B --- C , A-B --- C'
            v1 = coordC - coordB
            v2 = cC - coordB
            m = bRotation.get_rot_mat4v(v1, v2)
            mat.append(m)
            
        # do the transformation
        model = xyzModel()
        for m in mat:
            t = copy.deepcopy(xyz)
            t.transform(m, frg[1])
            model.extend(t)
        model.dump()    
        return model
        
        
    def mapping(self):
        pass
        
    def mk_frg(self, mol, sitelist):
        """ a list of mole """
            
        return
    def mk_frg(self, mol, frglist):
        """ split the mole in to fragments and group up"""
        
        return

    def regular_polygon(self, n=3, rad=3.0, theta=0.0, theta0=0.0, mode="circle"):
        """ build regular polygonal distribution """
        if n < 3:
            print "n < 3 in regular_polygon (??) !!!"
            exit(1)
        theta = theta / 180.0 * np.pi
        theta0 = theta0 / 180.0 * np.pi
        dt = 2.0 * np.pi / n
        dist = np.sqrt(2.0 * rad * rad * (1 - np.cos(dt)))
        print dist

        if mode == "polygon":
            dist = rad
            rad = dist / np.sqrt(2.0*(1-np.cos(dt)))

        bis_ang = np.arccos(dist/rad/2)
        
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
        for vec in edges:
            v = bRotation.rotate_ar_vector(vec, axis, theta+bis_ang)
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
        mt1 = bRotation.get_trans_mat4(-origin)

        for i in xrange(n_sides):
            # print vec, branches[i], origin, vertices[i]
            mt2 = bRotation.get_trans_mat4(vertices[i])

            mr = bRotation.get_rot_mat4v(vec, branches[i])
 
            m = np.dot(mr, mt1)

            m = np.dot(mt2, m)
            
            #m = bRotation.get_mat4t(vec, branches[i], origin, vertices[i]) 
            mat.append(m)
        self.polygon['mapping'] = mat
        return
        
    def print_polygon(self):
        fp = open("polygon.xyz", "w")
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
    
    def build(self, config):
    
        # mole
        core = config['mole']['core']
        direction = config['mole']['direction']
        
        # shape
        n_sides = config['polygon']['n_sides']
        t0 = config['polygon']['theta0']
        mode = config['polygon']['mode']
        # min max step
        rmin = config['radii']['min']
        rmax = config['radii']['max']
        rsize = config['radii']['size']
        tmin = config['theta']['min']
        tmax = config['theta']['max']
        tsize = config['theta']['size']
        # fies
        xyzfile = config['file']['xyz']
        #
        xyz = xyzSingle()
        xyz.read(filename=xyzfile)
        origin, vec = xyz.set_info(origin=core, direction=direction)
        polygon = bPolygon()
        
        rnum = int((rmax-rmin) / rsize)
        tnum = int((tmax-tmin) / tsize)  
        # 
        model = xyzModel()
        for i in xrange(rnum):
            for j in xrange(tnum):
                r = rmin + i * rsize
                t = tmin + j * tsize        
                polygon.regular_polygon(n=n_sides, rad=r, theta=t, theta0=t0, mode=mode)
                polygon.mapping(origin, vec)
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
    # template
    os.chdir("tmp")
    xyz = xyzSingle()
    xyz.read(filename="10.xyz")
    pes = bPESGrid()
    # pes.radical(xyz)
    pes.angle(xyz)
    # origin, vec = xyz.set_info(origin=0, direction=1)
    # polygon = bPolygon()
    # polygon.regular_polygon(n=5, rad=8.0, theta=3.0, theta0=0.0)
    # polygon.mapping(origin, vec)
    # newxyz = polygon.images(xyz)
    # newxyz.dump()
    
    # polygon.print_polygon()

    # line = input("enter working directory: \n > ")
    # os.chdir(line)

    # xyzfile = input("xyz file name: \n > ")
    
    # core, direction = input("enter X A number id: \n > ")

    
    # config = {}
    # config['radii'] = {'min': 3, 'max': 6, 'size': 0.1}
    # config['theta'] = {'min': -90, 'max': 90, 'size': 5}
    # config['polygon'] = {'n_sides': 3, 'theta0': 8, 'mode': 'polygon'}
    # config['mole'] = {'core': core, 'direction': direction}
    # config['file'] = {'xyz': xyzfile}
    # polygon = bPolygon()
    # polygon.build(config)
    
    
    
    
    # origin = np.array([0, 1, 1])
    # vec = np.array([3,3,3])
    # b = bPolygon()
    # b.regular_polygon(n=3, rad=3.0, theta=0.0, theta0=0.0)
    # b.mapping(origin, vec)

    
    
