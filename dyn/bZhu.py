com#! /usr/bin/env python


#
# the algorithm based on Zhu-Nakamura theory
#

class bZhu:
    def __init__(self):

        return

    def 

    def get_geom(self):
        """ read gaussian geom data """      
        filename = self.files['log']
        fp = open(filename, "r")
        n_atom = self.dim['n_atom']
        
        pat1 = re.compile("Input orientation:")
        pat2 = re.compile("Standard orientation:")
        pat3 = re.compile("Z-Matrix orientation:")
    
        std_pos = 0
        ipt_pos = 0
        z_pos = 0
        geom = [{} for i in xrange(n_atom)]
        while True:
            line = fp.readline()
            if line == "":
                break
            m = pat1.search(line)
            if m is not None:
                ipt_pos = fp.tell()
            m = pat2.search(line)
            if m is not None:
                std_pos = fp.tell()
            m = pat3.search(line)
            if m is not None:
                z_pos = fp.tell()
                
        # read geom
        if std_pos != 0:
            print "Read Standard orientation: Find Geometry"
            fp.seek(std_pos)
        elif ipt_pos != 0:
            print "Read Input orientation: Find Geometry"
            fp.seek(ipt_pos)
        elif z_pos != 0:
            print "Z-matrix orientation: Find Geometry.. Jujar"
            fp.seek(z_pos)
        else:
            print "FAILED TO READ GAUSSIAN FREQ. CHECK LOG FILE"
            exit(1)
        # start to read geom
        # jump four lines  
        for i in xrange(4):
            line = fp.readline()
        # read coord.
        for i in xrange(n_atom):
            coord = [0.0 for j in xrange(3)]
            line = fp.readline() 
            record = line.split()
            atom_number = int(record[1])
            coord[0] = float(record[3])
            coord[1] = float(record[4])
            coord[2] = float(record[5])
            atom = {'atom_number': atom_number, 'coord': coord}
            geom[i] = atom
        fp.close()

        # dump xyz file for check
        file_xyz = open("check.xyz", "w")
        n_atom = self.dim['n_atom']
        print >>file_xyz, "%10d\n" % n_atom
        for atom in geom:
            m = elements.elements()
            atom_label = m.number2label(atom['atom_number'])
            print >>file_xyz, "%10s%15.8f%15.8f%15.8f" % (atom_label,
                                                          atom['coord'][0],
                                                          atom['coord'][1],
                                                          atom['coord'][2])
        file_xyz.close()

        
        return        

# suppose direction is observed.

# find beta

# estimate minimum

# calc. alpha

# calc. PZN
