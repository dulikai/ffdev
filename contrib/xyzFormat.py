#! /usr/bin/python

# the xyz format from vmd, usually unsatisfactory
#
# xyz format: single..
#

# simplely read in a line, if HNXXX, then change it
#

import os
import copy

import shutil

class xyzFormat:
    config = {}
    config['mydir'] = 'ppp'
    
    table = {}
    table['list'] = ['h', 'c', 'n', 'o', 'f']
    table['set'] = set(table['list'])
    
    def __init__(self):
        self.model = []
        return
        
    def read_once(self, fpin):
        """ read in xyz format in angstrom unit """    
        # read number of atom
        line = fpin.readline().strip()  
        # check the end of file
        if line == "":
            return 0
        n_atom = int(line.split()[0])        
        title = fpin.readline().strip()
        geom = []
        atom_name = []
        # read a mole
        for i in xrange(n_atom):
            line = fpin.readline()
            rec = line.split()
            name, x, y, z= rec[0:4]            
            coord = [float(x),float(y),float(z)]
            atom_name.append(name)
            geom.append(coord)
        mol = {'title':title, 'n_atom': n_atom, 'geom':geom, 'name': atom_name}
        return mol        
    
    def read(self, filename="dump.xyz"):
        """ read in a list of xyz """
        fp = open(filename, "r")
        while True:
            mol = self.read_once(fp)
            if mol == 0:
                break
            self.model.append(mol)
        fp.close()    
        return
        
                
    def modify_once(self, mol):
        """ atom name etc. """
        n_atom = mol['n_atom']
        name = mol['name']
        myname = ["O" for i in xrange(n_atom)]
        mydict = {}
        # print self.table['set']
        for i in xrange(n_atom):
            atom_name = name[i]
            t0 = atom_name in mydict.keys()
            t1 = atom_name[0:1].lower() in self.table['set']
            t2 = atom_name[0].lower() in self.table['set']
            if t0:
                myname[i] = mydict[atom_name]
            elif t1:
                myname[i] = atom_name[0:1].lower()
                mydict[atom_name] = myname[i]
            elif t2:
                myname[i] = atom_name[0].lower()
                mydict[atom_name] = myname[i]
            elif t1 and t2:
                print "conflict element name for type %s" % atom_name
                line = raw_input("enter element name manually:\n > ")
                myname[i] = line.strip().lower()
                mydict[atom_name] = myname[i]                
            else:
                print "impossible element name for type %s" % atom_name
                line = raw_input("enter element name manually:\n > ")
                myname[i] = line.strip().lower()
                mydict[atom_name] = myname[i]                
            # print "%s" % myname[i]            
        mol['name'] = myname    
        return
        
    def modify(self):
        """ mod all """
        model = self.model
        for mol in model:
            self.modify_once(mol)
        return
        
    def dump_once(self, mol, filename="punch.xyz"):
        """ write xyz format in angstrom unit """    
        if os.path.isfile(filename):
            print "OVERWRITE: %s" % filename
        fp = open(filename, "w")
        title = mol['title']
        n_atom = mol['n_atom']
        geom = mol['geom']
        name = mol['name']
        print >>fp, "%10d" % (n_atom)
        print >>fp, title
        for i in xrange(n_atom):
            pos = geom[i]
            atom_name = name[i]
            print >>fp, "%10s%12.6f%12.6f%12.6f" % (atom_name, pos[0], pos[1], pos[2])
        fp.close()    
        
        return
    
    def dump_many(self, filename="punch.xyz"):
        """ write all xyz files """
        mydir = self.config['mydir']
        if os.path.isdir(mydir):
            shutil.rmtree(mydir)
        os.mkdir(mydir)
        os.chdir(mydir)
        mystr = filename.split(".")
        prefix = mystr[0]
        model = self.model
        for i in xrange(len(model)):
            mol = model[i]
            myfile = filename + "." + str(i)
            self.dump_once(mol, myfile)
            print "MOL %d done" % i
        os.chdir("../")    
        return    
 
    def print_mol(self, mol):
        """ print mol into string """
        mystr = ""
        title = mol['title']
        n_atom = mol['n_atom']
        geom = mol['geom']
        name = mol['name']
        mystr += "%10d\n" % n_atom
        mystr += title + "\n"
        for i in xrange(n_atom):
            pos = geom[i]
            atom_name = name[i]
            mystr += "%10s%12.6f%12.6f%12.6f\n" % (atom_name, pos[0], pos[1], pos[2])
        return mystr
        
    def dump(self, filename="punch.xyz"):
        """ write all xyz files """
        fp = open(filename, "w")
        model = self.model
        for mol in model:
            mystr = self.print_mol(mol)
            print >>fp, "%s" % mystr,
        fp.close()    
        return    
 
if __name__ == "__main__":
    
    xyz = xyzFormat()
    xyz.read(filename="./dump.xyz")
    xyz.modify()
    xyz.dump()
    xyz.dump_many()




