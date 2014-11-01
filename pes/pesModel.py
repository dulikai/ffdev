#! /usr/bin/python

# here define one of the most important format 
# xyz format: single..

import os
import copy
import numpy as np


def read_pes( filename="pes.dat"):
    """ read in one set of pes """
    au2kcal = 627.5
    pes = {}
    fp = open(filename, "r")
    n_site = int(fp.readline())
    title = fp.readline()
    sites = []
    for i in xrange(n_site):
        line = fp.readline()
        val = [float(v) for v in line.split()] 
        s1 = (val[0] - val[1] - val[2]) * au2kcal
        s2 = (val[0] - val[3] - val[4]) * au2kcal
        body = [s1, s2]
        sites.append(body)
    fp.close()
    pes['sites'] = sites
    pes['n_site'] = n_site
    pes['title'] = title
        
    fp.close()
    return pes

def read_vars(filename="vars.dat"):
    fp = open(filename, "r")
    vars = []
    while True:
        line = fp.readline().strip()
        if line == "":
            break
        vars.append(line)
    return vars

def dump_pes(pes, vars):
    n_site = pes['n_site']
    sites = pes['sites']
    title = pes['title']
    fp = open("dump.dat", "w")
    print >>fp, n_site
    print >>fp, title
    for i in xrange(n_site):
        print >>fp, "%s%12.6f%12.6f" % (vars[i], sites[i][0], sites[i][1])
    fp.close()
    return
if __name__ == "__main__":
    pes = read_pes(filename="pes.dat")
    vars = read_vars(filename="vars.dat")
    dump_pes(pes, vars)

    
       

        
        
# if __name__ == "__main__":
    # pes = pesSingle()
    # pes.read()
    # pes.dump()


