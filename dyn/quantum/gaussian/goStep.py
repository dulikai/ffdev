#! /usr/bin/env python

import os
import subprocess
import shutil

from SetupEnv import *
from CreateInp import *
from ffRunner import *
from Parser import *

class goStep():
    """
    singlet state dynamics
    """
    def __init__(self, config = {}):
        """
        parameter set.
        """
        self.config = config
        
        return

        
    def setup(self):
        """
        check interface, make directory,
        setup the job exec. environment
        """
        s = SetupEnv(self.config)
        s.setup()
        return
        
    def create(self):
        """ generate input file """
        c = CreateInp(self.config)
        c.modify(jobtype = "dft")
        c.wrt_input()
        return
        
    def runner(self):
        """ exec qc calc. """
        ff = ffRunner(self.config)
        ff.execjob()
        return
        
    def parser(self):
        """ parser ff calculation output """
        p = Parser()
        p.get()
        return
        
    def collect(self):
        """ collect output info. """
        
        return
        
        
    def step(self):
        """
        build directory, generate input, exec job, parser info, collect data
        """
        self.setup()
        
        self.create()
        
        self.runner()
        
        self.parser()
        
        self.collect()
        
        return
        
        
    def collect(self):
        """
        simply clean up the tmp dat. and so on.
        """
        #   Go back to directory of dynamics work
        #   Copy results of QM calculations
        
        #   Go back to directory of dynamics work
        os.chdir(self.directory['root'])
                            
        #   Copy results of QM calculations 
        sourcePath = self.directory['work']
        sourceFile = sourcePath + '/' + 'qm_results.dat'
        destPath = './'
        shutil.copy2(sourceFile, destPath)

        sourceFile = sourcePath + '/' + 'qm_other.dat'
        destPath = './'
        shutil.copy2(sourceFile, destPath)

        print 'Finish QC calculation'          
        
        return
        

if __name__ == "__main__":
    config = tools.load_data("config.in")
    q = goStep(config)
    q.step()
    
    
        
        

