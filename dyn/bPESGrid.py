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
    

    def ang1(self, mol):
         # suppose the xyz geom is given..
        
        frg = [[0,1,2,3,4], [5]]
        ndx = [[0], [4], [5]]        
        tbl = {"start": 90.0, "end": 270.0, "size": 1.0}
        
        # build ndx center for angle moving..
        cA = xyz.get_center(ndx[0])
        cB = xyz.get_center(ndx[1])
        cC = xyz.get_center(ndx[2])
        # mapping   
        # setup A B C at the origin
        # C --- B --- A --> x axis
        # frag 1 no rotation
        rAB = np.linalg.norm(cB-cA)
        oA = np.array([-rAB, 0.0, 0.0])
        oB = np.zeros(3)
        # then we can build a series of A, B, C pairs,
        # for which, the radical scan is done..
        start = tbl['start'] / 180.0 * np.pi
        end = tbl['end'] / 180.0 * np.pi
        size = tbl['size'] / 180.0 * np.pi
        n_points = int((end-start)/size)
        pairs = []
        fp = open("test.xyz", "w")
        for i in xrange(n_points):
            # suppose in xy plane
            theta = start + size * i
            print i, theta / np.pi * 180.0, oC
            rBC = np.linalg.norm(cC-cB)
            x = rBC * np.cos(np.pi-theta)
            y = rBC * np.sin(np.pi-theta)
            oC = np.array([x, y, 0.0])
            pairs.append([oA, oB, oC])
            print >>fp, "3"
            print >>fp, ""
            print >>fp, "O %12.6f%12.6f%12.6f" % (oA[0], oA[1], oA[2])
            print >>fp, "C %12.6f%12.6f%12.6f" % (oB[0], oB[1], oB[2])
            print >>fp, "O %12.6f%12.6f%12.6f" % (oC[0], oC[1], oC[2])
            # print oA, oB, oC
        fp.close()
        return
        
    def angle(self, mol):
        """
        frag 1 is fixed by default..
        """
         # suppose the xyz geom is given..
        
        frg = [[0,1,2,3,4], [5]]
        ndx = [[0], [4], [5]]        
        tbl = {"start": 90.0, "end": 270.0, "size": 1.0}
        
        # build ndx center for angle moving..
        cA = xyz.get_center(ndx[0])
        cB = xyz.get_center(ndx[1])
        cC = xyz.get_center(ndx[2])
        
        # mapping   
        # setup A B C at the origin
        # C --- B --- A --> x axis
        # frag 1 no rotation
        rAB = np.linalg.norm(cB-cA)
        oA = np.array([-rAB, 0.0, 0.0])
        oB = np.zeros(3)
        
        # then we can build a series of A, B, C pairs,
        # for which, the radical scan is done..
        start = tbl['start'] / 180.0 * np.pi
        end = tbl['end'] / 180.0 * np.pi
        size = tbl['size'] / 180.0 * np.pi
        n_points = int((end-start)/size)
        pairs = []
        for i in xrange(n_points):
            # suppose in xy plane
            theta = start + size * i
            rBC = np.linalg.norm(cC-cB)
            x = rBC * np.cos(np.pi-theta)
            y = rBC * np.sin(np.pi-theta)
            oC = np.array([x, y, 0.0])
            pairs.append([oA, oB, oC])
            # print oA, oB, oC
        # build transformation matrix @
        mat = []
        for i in xrange(n_points):
            oA = pairs[i][0]
            oB = pairs[i][1]
            oC = pairs[i][2]
            # A-B --- C , A'-B' --- C'
            oBA = oA - oB
            oBC = oC - oB
            cBA = cA - cB
            cBC = cC - cB
            # build matrix
            # frag 0
            mAt = bRotation.get_trans_mat4(oB-cB)
            mAr = bRotation.get_rot_mat4v(cBA, oBA)
            # frag 1
            mBt = bRotation.get_trans_mat4(oC-cC)
            mBr = bRotation.get_rot_mat4v(cBC, oBC)
            # summary
            mat.append([mAt, mAr, mBt, mBr])
        # do the transformation
        model = xyzModel()
        for m in mat:
            t = copy.deepcopy(xyz)
            excl = ndx[1]  # B site
            t.transform(m[0], frg[0])  # trans.
            t.transform(m[1], frg[0], excl) # rot.
            t.transform(m[2], frg[1])  # trans.
            excl = ndx[2]  # C site
            t.transform(m[3], frg[1], excl)  # rot.
            model.extend(t)
        model.dump()    
        return model
        
        
    def dihedral1(self, mol):
        """ 
        the dihedral phi, is A-X --- Y-M dihedral.
        that is X---Y is at x-axis.

        the dihedral phi, is A-X --- Y-M dihedral.
        In this transformation, we maintain the theta1 theta2, and the distance between X/Y
        SO the dihedral is possible to solve exactly.
        THIS is: move the atom A/M(refer as Z) to definitive position.
        DEFINE Z(x,y,z)
        IN our coordinate system, YXA or XYM. 
        TO satisfy the above constraint, x coord is known to be fixed; y/z is movable.
        THEY(x,y) satisfy: z^2+y^2 = x^2; z/y = tan[theta()] KNOWN x is constant.
		"""
        frg = [[0,1,2,3,4], [5,6,7,8,9]]
        ndx = [[0], [4], [9], [5]]        
        tbl = {"start": 90.0, "end": 270.0, "size": 1.0}
        
        # build ndx center for DIHEDRAL moving..
        cA = xyz.get_center(ndx[0])
        cX = xyz.get_center(ndx[1])
        cY = xyz.get_center(ndx[2])
        cM = xyz.get_center(ndx[3])
 
        # mapping   
        # Y at the origin, and X --- Y at x axis.
        # M--Y--X--A
        # Y center
        oY = np.zeros(3)
        # X center
        vYX = cX - cY
        rXY = np.linalg.norm(vYX)
        oX = np.array([rXY, 0.0, 0.0])
        
        # theta1 is A-X-Y angle
        vXY = cY - cX; vXA = cA - cX;
        rXY = np.linalg.norm(vXY)
        rXA = np.linalg.norm(vXA)
        costheta1 = np.dot(vXY, vXA) / (rXY * rXA)
        # theta2 is X-Y-M angle
        vYX = -vXY; vYM = cM - cY;
        rYX = rXY
        rYM = np.linalg.norm(vYM)
        costheta2 = np.dot(vYX, vYM) / (rYX * rYM)
        # theta 1 2
        theta1 = np.arccos(costheta1)
        theta2 = np.arccos(costheta2)
        
        # A center
        z = 0.0
        y = -rXA * np.sin(theta1)
        x = rXY - rXA * np.cos(theta1)
        oA = np.array([x, y, z])
        
        start = tbl['start'] / 180.0 * np.pi
        end = tbl['end'] / 180.0 * np.pi
        size = tbl['size'] / 180.0 * np.pi
        n_points = int((end-start)/size)
        pairs = []
        fp = open("test.xyz", "w")
        for i in xrange(n_points):
            # suppose in Y, X point in x axis
            phi = start + size * i
            # M center
            x = rYM * np.cos(theta2)
            y = x * np.cos(phi)
            z = x * np.sin(phi)
            oM = np.array([x, y, z])
            pairs.append([oA, oX, oY, oM])
            print >>fp, "4"
            print >>fp, ""
            print >>fp, "H %12.6f%12.6f%12.6f" % (oA[0], oA[1], oA[2])
            print >>fp, "O %12.6f%12.6f%12.6f" % (oX[0], oX[1], oX[2])
            print >>fp, "O %12.6f%12.6f%12.6f" % (oY[0], oY[1], oY[2])
            print >>fp, "H %12.6f%12.6f%12.6f" % (oM[0], oM[1], oM[2])
        fp.close()    
        return
        
        
    def dihedral(self, mol):
        """ 
        the dihedral phi, is A-X --- Y-M dihedral.
        that is X---Y is at x-axis.

        the dihedral phi, is A-X --- Y-M dihedral.
        In this transformation, we maintain the theta1 theta2, and the distance between X/Y
        SO the dihedral is possible to solve exactly.
        THIS is: move the atom A/M(refer as Z) to definitive position.
        DEFINE Z(x,y,z)
        IN our coordinate system, YXA or XYM. 
        TO satisfy the above constraint, x coord is known to be fixed; y/z is movable.
        THEY(x,y) satisfy: z^2+y^2 = x^2; z/y = tan[theta()] KNOWN x is constant.
		"""
        frg = [[0,1,2,3,4], [5,6,7,8,9]]
        ndx = [[0], [4], [9], [5]]        
        tbl = {"start": -90.0, "end": 90.0, "size": 1.0}
        
        # build ndx center for DIHEDRAL moving..
        cA = xyz.get_center(ndx[0])
        cX = xyz.get_center(ndx[1])
        cY = xyz.get_center(ndx[2])
        cM = xyz.get_center(ndx[3])
 
        # mapping   
        # Y at the origin, and X --- Y at x axis.
        # M--Y--X--A
        # Y center
        oY = np.zeros(3)
        # X center
        vYX = cX - cY
        rXY = np.linalg.norm(vYX)
        oX = np.array([rXY, 0.0, 0.0])
        
        # theta1 is A-X-Y angle
        vXY = cY - cX; vXA = cA - cX;
        rXY = np.linalg.norm(vXY)
        rXA = np.linalg.norm(vXA)
        costheta1 = np.dot(vXY, vXA) / (rXY * rXA)
        # theta2 is X-Y-M angle
        vYX = -vXY; vYM = cM - cY;
        rYX = rXY
        rYM = np.linalg.norm(vYM)
        costheta2 = np.dot(vYX, vYM) / (rYX * rYM)
        # theta 1 2
        theta1 = np.arccos(costheta1)
        theta2 = np.arccos(costheta2)
        
        # A center
        z = 0.0
        y = -rXA * np.sin(theta1)
        x = rXY - rXA * np.cos(theta1)
        oA = np.array([x, y, z])
        print np.sqrt(x*x+y*y+z*z)
        print np.linalg.norm(oA-oX)
        start = tbl['start'] / 180.0 * np.pi
        end = tbl['end'] / 180.0 * np.pi
        size = tbl['size'] / 180.0 * np.pi
        n_points = int((end-start)/size)
        pairs = []
        fp = open("test.xyz", "w")
        for i in xrange(n_points):
            # suppose in Y, X point in x axis
            phi = start + size * i
            # M center
            x = rYM * np.cos(theta2)
            rad = rYM * np.sin(theta2)
            y = rad * np.cos(phi)
            z = rad * np.sin(phi)
            oM = np.array([x, y, z])
            pairs.append([oA, oX, oY, oM])
            print >>fp, "4"
            print >>fp, ""
            print >>fp, "H %12.6f%12.6f%12.6f" % (oA[0], oA[1], oA[2])
            print >>fp, "O %12.6f%12.6f%12.6f" % (oX[0], oX[1], oX[2])
            print >>fp, "O %12.6f%12.6f%12.6f" % (oY[0], oY[1], oY[2])
            print >>fp, "H %12.6f%12.6f%12.6f" % (oM[0], oM[1], oM[2])
        fp.close()    
        # build transformation matrix @
        mat = []
        for i in xrange(n_points):
            oA = pairs[i][0]
            oX = pairs[i][1]
            oY = pairs[i][2]
            oM = pairs[i][3]
            # A-X Y-M to oXA oYM
            oXA = oA - oX
            cXA = cA - cX
            oYM = oM - oY
            cYM = cM - cY
            # build matrix
            # frag 0
            mAt = bRotation.get_trans_mat4(oX-cX) # X center
            mAr = bRotation.get_rot_mat4v(cXA, oXA)
            # frag 1
            mBt = bRotation.get_trans_mat4(oY-cY) # Y center
            mBr = bRotation.get_rot_mat4v(cYM, oYM)
            print np.linalg.norm(oA)

            # summary
            mat.append([mAt, mAr, mBt, mBr])
        # do the transformation
        model = xyzModel()
        for m in mat:
            t = copy.deepcopy(xyz)
            excl = ndx[1]  # X site
            t.transform(m[0], frg[0])  # trans.
            t.transform(m[1], frg[0], excl) # rot.
            t.transform(m[2], frg[1])  # trans.
            excl = ndx[2]  # Y site
            t.transform(m[3], frg[1], excl)  # rot.
            model.extend(t)
        model.dump()    
        return
        
        
	def __gen_phi_geom_mol(self, grid):
		""" move one molecule coordinate, ie. one frame """
		#work_mol
		mol = self.work_mol
		natom = mol['natom']
		atom = mol['atom']	
		title = mol['title']
		point = grid['coord']
		index = grid['index']
		title = "%s %s" % (title, index)
		cX = self.keypoint[1]
		cY = self.keypoint[2]
		oldcoord = self.keypoint[3]
		if cX[0] - cY[0] < 0.0:
			oldcoord = self.keypoint[0]
		xaxis = np.array([oldcoord[0], 0.0, 0.0])
		oldcoord = np.subtract(oldcoord, xaxis)
		point = np.subtract(point, xaxis)
		theta, vdir = self.__get_axis_angle(oldcoord, point)
		new_atom = []		
		cmin = 1.0e-4
		for i in range(natom):
			record = atom[i]
			if record['iflag'] != 0 and theta > cmin:
				# A-X --- Y-M, the M and X vector
				coord = self.point_rotate_ar_vector(record['coord'], vdir, theta)
			else:
				coord = record['coord']
			tmp_record = {'name': record['name'], 'coord': coord, 'iflag': record['iflag'], 'index': grid['index'] }
			new_atom.append(tmp_record)
		self.tmp_mol = {'natom':natom, 'atom':new_atom, 'title':title}	
		return	        
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
    xyz.read(filename="ch3br2.xyz")
    pes = bPESGrid()
    # pes.radical(xyz)
    # pes.angle(xyz)
    pes.dihedral(xyz)
    # pes.ang1(xyz)
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

    
    
