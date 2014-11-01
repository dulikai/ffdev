#! /usr/bin/python

# here define one of the most important format 
# xyz format: single..

import os
import copy
import numpy as np

class XYZGeom:
    
    """Hold a cartesian geometry.

    This class provides methods to calculate useful structural quantities
    such as bond lengths, bond angles, dihedral angles.

    The geometry can be written out in XYZ or MNDO formats.

    """

    def __init__(self):
        self.noOfAtoms = 0
        self.atomicNumbers = []
        self.atomCoords = []

    def addAtom(self, atomicNumber, x, y, z):
        """Adds an atom and cartesian coordinates to the geometry."""
        self.noOfAtoms += 1
        self.atomicNumbers.append(int(atomicNumber))
        newCoord = XYZPoint(float(x), float(y), float(z))
        self.atomCoords.append(newCoord)

    def getBondLength(self, a1, a2):
        """Return the distance between atom numbers a1 and a2.

        Atoms are numbered from zero.

        """
        diff = self.atomCoords[a1] - self.atomCoords[a2]
        return diff.getNorm()

    def getBondAngle(self, a1, a2, a3, degrees=True):
        """Return the bond angle a1-a2-a3.

        The angle is defined by the vectors a1-a2 and a2-a3.
        Atoms are numbered from zero.

        degrees - if true, return angle in degrees, else radians.

        """
        # Based on the bond angle routine from geoman.f90 by Eduardo Fabiano
        # vector 1 = 1->2
        v1 = self.atomCoords[a2] - self.atomCoords[a1]
        v1 = v1.normalise()

        # vector 2 = 2->3
        v2 = self.atomCoords[a3] - self.atomCoords[a2]
        v2 = v2.normalise()

        # dot product
        dotp = v1.dotProduct(v2)
 
        # angle in radians
        ang = math.pi - math.acos(dotp)
        if degrees:
            ang *= (180.0 / math.pi)
        return ang

    def getOutOfPlaneAngle(self, a1, a2, a3, a4, degrees=True):
        """Return the out of plane angle a1-a2(-a4)-a3.

        The angle is defined between the plane a1-a2-a3 and the vector a2-a4.
        Atoms are numbered from zero.

        degrees - if true, return angle in degrees, else radians.

        """
        # vector 2->1
        v21 = self.atomCoords[a1] - self.atomCoords[a2]
        v21 = v21.normalise()

        # vector 2->3
        v23 = self.atomCoords[a3] - self.atomCoords[a2]
        v23 = v23.normalise()

        # vector 2->4
        v24 = self.atomCoords[a2] - self.atomCoords[a4]
        v24 = v24.normalise()

        # vector product v21^v23
        w = v21.crossProduct(v23)
        w = w.normalise()

        # dot product
        dotp = w.dotProduct(v24)

        # angle in radians
        # using sin as we are actually calculating the angle
        # to the normal of the plane rather than the plane itself.
        ang = math.asin(dotp)        

        if degrees:
            ang *= (180.0 / math.pi)
        return ang   
                   
    def getDihedralAngle(self, a1, a2, a3, a4, degrees=True, signed=True):
        """Return the dihedral angle a1-a2-a3-a4.

        The angle is defined between the planes a1-a2-a3 and a2-a3-a4.
        Atoms are numbered from zero.

        degrees - if true, return angle in degrees, else radians.

        signed - if true, return a signed dihedral angle according to
                 MNDO conventions

        """
        # Based on the dihedral routine from geoman.f90 by Eduardo Fabiano
        # vector 1 = 1->2
        v1 = self.atomCoords[a2] - self.atomCoords[a1]
        v1 = v1.normalise()

        # vector 2 = 2->3
        v2 = self.atomCoords[a3] - self.atomCoords[a2]
        v2 = v2.normalise()

        # vector 3 = 3->4
        v3 = self.atomCoords[a4] - self.atomCoords[a3]
        v3 = v3.normalise()

        # vector product 1 = v1^v2
        w1 = v1.crossProduct(v2)
        w1 = w1.normalise()

        # vector product 2 = v3^v2
        w2 = v3.crossProduct(v2)
        w2 = w2.normalise()

        # dot product
        dotp = w1.dotProduct(w2)

        # angle in radians
        ang = math.pi - math.acos(dotp)
        if signed:
            # the dot product between the 1->2 bond and the normal to
            # the plane 2-3-4 is used to determine whether the dihedral
            # angle is clockwise or anti-clockwise
            # if the dot product is positive, the rotation is clockwise
            #
            # Taken from ChemShell FRAG_dihedral in interface.c
            #
            rotdir = v1.dotProduct(w2)
            if rotdir > 0:
                ang = -ang
        if degrees:
            ang *= (180.0 / math.pi)
        return ang        


    def getSphereRadius(self):
        """Return the radius of molecules and et al..

        Atoms are numbered from zero.

        """
        noOfAtoms = self.noOfAtoms
        averCoord = XYZPoint(0.0, 0.0, 0.0)
        for i in xrange(noOfAtoms):
            averCoord = averCoord + self.atomCoords[i]
        averCoord = averCoord.scale(1.0/noOfAtoms)

        radius = 0.0
        for i in xrange(noOfAtoms):
            diff = self.atomCoords[i] - averCoord
            radius += diff.getNorm()
        radius /= noOfAtoms 

        return radius


    def writeXYZGeom(self, title='', symbols=True):
        """Return the geometry in XYZ format

        title: sets title line
        symbols: if True, write atomic symbols instead of numbers

        Return a list containing the XYZ formatted geometry

        """
        if symbols:
            p = PeriodicTable()
        xyzGeom = []
        xyzGeom.append("%5s" % self.noOfAtoms)
        xyzGeom.append(title)
        for i in range(self.noOfAtoms):
            if symbols:
                atom = p.getSymbol(self.atomicNumbers[i])
            else:
                atom = self.atomicNumbers[i]
            geomLine = "%3s %20.10f %20.10f %20.10f" % (atom,
                                                        self.atomCoords[i].x,
                                                        self.atomCoords[i].y,
                                                        self.atomCoords[i].z)
            xyzGeom.append(geomLine)
        return xyzGeom



class PeriodicTable:
    
    """Hold a list of element names and a dictionary of atomic numbers."""

    def __init__(self):
        self.symbols = ["-", "H", "He",
                        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                        "K", "Ca",
                        "Sc", "Ti", "V", "Cr", "Mn",
                        "Fe", "Co", "Ni", "Cu", "Zn",
                        "Ga", "Ge", "As", "Se", "Br", "Kr",
                        "Rb", "Sr",
                        "Y", "Zr", "Nb", "Mo", "Tc",
                        "Ru", "Rh", "Pd", "Ag", "Cd",
                        "In", "Sn", "Sb", "Te", "I", "Xe",
                        "Cs", "Ba",
                        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                        "Hf", "Ta", "W", "Re",
                        "Os", "Ir", "Pt", "Au", "Hg",
                        "Tl", "Pb", "Bi", "Po", "At", "Rn",
                        "Fr", "Ra",
                        "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
                        "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                        "Rf", "Db", "Sg", "Bh",
                        "Hs", "Mt", "Ds", "Rg"]
        self.atomicNumbers = {}
        for i in range(1, len(self.symbols)):
            self.atomicNumbers[self.symbols[i].lower()] = i

    def getSymbol(self, atomicNumber):
        """Return an atomic symbol given an atomic number"""
        if atomicNumber < 1 or atomicNumber >= len(self.symbols):
            return None
        else:
            return self.symbols[atomicNumber]

    def getZ(self, symbol):
        """Return an atomic number given an atomic symbol"""
        mysymbol = symbol.lower()
        if self.atomicNumbers.has_key(mysymbol):
            return self.atomicNumbers[mysymbol]
        else:
            return None


            

class xyzPoint:
    """Hold a 3 dimensional Cartesian coordinate vector."""
    name = "Ar"
    pos = np.zeros(3)

    def __init__(self, coord=np.zeros(3), name="C"):
        """Initialise coordinates."""
        self.pos = coord[0:3]
        self.name = name

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
            
            
class xyzSingle:
    """ a single xyz geom. """
    
    def __init__(self):
        self.sites = []
        self.n_site = 0
        self.center_id = 0
        self.dir_id = [0, 1]
        self.title = ""
        return
        
    def __dir__(self):
        return ["sites", "n_site"]

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
        sites = self.give(keyword="sites")
        for ibody in sites:
            ibody.pos = np.zeros(3)
            ibody.name = "C"
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
        
    def read(self, filename="efp.xyz"):
        """
        read cart. coordinate in xyz format
        """
        # xyz file
        fp = open(filename, "r")
        n_site = int(fp.readline())
        self.title = fp.readline()
        sites = [xyzPoint() for i in xrange(n_site)]
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
        # done
        return
            
        
    def extend(self, xyz):
        """ append another xyz file """
        sites = xyz.give(keyword="sites")
        self.sites.extend(sites)        
        self.n_site += xyz.n_site
        return
        
    def bond_vector(self, i, j):
        """
        calculate bond vector i ---> j
        i,j start from 0
        """
        n_site = self.n_site    
        if i > n_site or i < 0:
            print "index error"
        if j > n_site or j < 0:
            print "index error"            
        pi = self.sites[i]; pj = self.sites[j]    
        vec = pj - pi
        return vec.pos
        
    def point(self, i):
        """
        return point i
        """
        n_site = self.n_site    
        if i > n_site or i < 0:
            print "index error"            
        point = self.sites[i];
        return point.pos
        
    def set_info(self, origin=0, direction=1):
        """ set the rotation scheme """
        self.center_id = origin
        self.dir_id = [origin, direction]
        p = self.point(origin)
        v = self.bond_vector(origin, direction)
        return p, v
        
    def transform(self, m):
        """ rigid body motion of this mole. """
        n_site = self.n_site
        for i in xrange(n_site):
            m4 = copy.deepcopy(m)
            # if i == self.center_id:
                # m4[0:3,0:3] = np.eye(3)
            print m4[3][0:3]    
            m4[3][0:3] += self.point(self.center_id)
            self.sites[i] = self.sites[i].transfrom(m4)
            # print m4
        return 
            
    def read_xyz(self):
        """ read in xyz format in angstrom unit """    
        filename = self.files['xyz']        
        fpin = open(filename, "r")
        # read number of atom
        line = fpin.readline()  
        n_atom = int(line.split()[0])        
        title = fpin.readline().strip()
        geom = []
        atom_name = []
        # read a mole
        for i in xrange(n_atom):
            line = fpin.readline()
            rec = line.split()
            name, x, y, z= rec[0:4]            
            pos = np.array([float(x),float(y),float(z)])
            point = xyzPoint(pos)
            atom_name.append(name)
            geom.append(coord)
        mol = {'title':title, 'n_atom': n_atom, 'geom':geom, 'name': atom_name}
        fpin.close()
        self.mol = mol
        return mol

    def dump(self, filename="dump.xyz"):
        """ write xyz format in angstrom unit """    
        if os.path.isfile(filename):
            print "OVERWRITE: %s" % filename
        fp = open(filename, "w")
        sites = self.sites
        n_site = self.n_site
        print >>fp, "%10d" % (n_site)
        print >>fp, ""
        for ibody in sites:
            pos = ibody.pos
            name = ibody.name
            print >>fp, "%10s%12.6f%12.6f%12.6f" % (name, pos[0], pos[1], pos[2])
        fp.close()    
        return 

        
    def __str__(self):
        """Return simple string representation of the xyz geom"""
        sites = self.sites
        n_site = self.n_site
        title = self.title
        mystr  = "%10d\n" % n_site
        mystr += "%s\n" % title
        for ibody in sites:
            pos = ibody.pos
            name = ibody.name
            mystr += "%10s%12.6f%12.6f%12.6f\n" % (name, pos[0], pos[1], pos[2])
        # print n_site, len(sites) 
        return mystr
        
        

class xyzModel():
    
    """Hold an XYZ file containing one or more XYZ structures.
    """
    model = []
    n_xyz = 0
    
    def __init__(self):
    
        return
    
    def load(self, filename=""):
    
        return
        
    def dump(self, filename="dump.xyz"):
        """ dump many xyz geom """
        fp = open(filename, "w")
        for geom in self.model:
            print >>fp, geom,
        fp.close()    
        return
        
    def extend(self, xyz):
        self.model.append(xyz)
        self.n_xyz += 1
        return
        
    def read(self, fileName):
        """Load in the entire file as a list."""
        self.fileName = fileName
        try:
            f=open(fileName, 'r')
            self.fileContents = f.readlines()
            f.close()
        except IOError:
            sys.stderr.write("Error: could not read in %s\n" % fileName)
            sys.exit(1)
        # Load in periodic table for future use
        self.pt = PeriodicTable()
        # Check if the file is a Molden file, and if so, remove everything
        # before the XYZ information
        # Some useful information and very basic checks
        firstAtomLine = self.fileContents[0].split()
        try:
            self.noOfAtoms = int(firstAtomLine[0])
        except ValueError:
            sys.stderr.write("Error: could not read no. of atoms from first "
                             + "line of XYZ file.\n")
            sys.exit(1)
        if self.noOfAtoms <= 0:
            sys.stderr.write("Error: claimed no. of atoms on first line of " +
                             "XYZ file is <= 0.\n")
        # Cursory check of contents
        if len(self.fileContents) % (self.noOfAtoms + 2) == 0:
            self.noOfGeoms = len(self.fileContents) / (self.noOfAtoms + 2)
        else:
            sys.stderr.write("Error: XYZ file length is not an integer " +
                             "multiple of claimed individual size.\n")
            sys.exit(1)

    def getGeom(self, geomNumber):
        """Return an XYZGeom object containing a Cartesian geometry.

        geomNumber - specifies which geometry in the file should be parsed
                     (numbered from 0)
        
        Return None if the XYZ data could not be parsed.

        """
        xyzGeom = XYZGeom()
        # We calculate first line assuming that geomNumber is numbered from 0
        lineNo = geomNumber * (self.noOfAtoms + 2)
        claimedAtoms = int(self.fileContents[lineNo].split()[0])
        if claimedAtoms != self.noOfAtoms:
            sys.stderr.write("Error: claimed no. of atoms is inconsistent.\n")
            return None
        # Ignore the title line
        lineNo += 2
        for i in range(self.noOfAtoms):
            XYZLine = self.fileContents[lineNo].split()
            try:
                an = int(XYZLine[0])
            except ValueError:
                an = self.pt.getZ(XYZLine[0])
                if an is None:
                    sys.stderr.write("Error: first xyz column should be an "
                                     "atomic number or symbol.\n")
                    return None
            try:
                x = float(XYZLine[1])
                y = float(XYZLine[2])
                z = float(XYZLine[3])
            except ValueError:
                sys.stderr.write("Error: columns 2-4 should contain xyz " +
                                 "data separated by spaces.\n")
                return None
            xyzGeom.addAtom(an, x, y, z)
            lineNo += 1
        return xyzGeom

    def writeSnapshots(self, snapshots, fileName):
        """Write a new XYZ file containing snapshots from the original file.

        snapshots - list of geometries to be written (numbered from 0)
        fileName - name of the output file containing the snapshots

        The file is copied as-is, with no further checking beyond that when
        the object was initialised.

        Return -1 if unsuccessful

        Usage note: snapshots at regular intervals can be taken using 'range'
          i.e.. snapshots = range([1st snapshot], [total geoms], [interval])

        """
        try:
            f=open(fileName, 'w')
        except IOError:
            sys.stderr.write("Error: could not open %s" +
                             "for writing.\n" % fileName)
            return -1
        for snap in snapshots:
            if snap > self.noOfGeoms - 1:
                sys.stderr.write("Warning: requested snapshot %i is too "
                                 "large - ignored.\n" % snap)
                continue
            startLine = snap * (self.noOfAtoms + 2)
            for i in range(self.noOfAtoms + 2):
                f.write(self.fileContents[startLine + i])
        f.close()

        
        
        
if __name__ == "__main__":
    xyz = xyzSingle()
    xyz.read()


