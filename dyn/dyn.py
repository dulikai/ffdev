#! /usr/bin/env python

import os
import math
import numpy as np

# this is a very very very simple md eigne 
# first @SDU by dulikai
# then  @qibebt dulikai
#
#
# aim at propagation of the nuclear position
#


struct body {
	mat_t rotmat;
	vec_t pos;
	vec_t vel;
	vec_t vel_old;
	vec_t angmom;
	vec_t angmom_old;
	vec_t force;
	vec_t torque;
	vec_t inertia;
	vec_t inertia_inv;
	double mass;
};

class sim_site:
    def __init__(self):
        zero3d = [0.0, 0.0, 0.0]
        self.pos = np.array(zero3d)
        self.pos_old = np.array(zero3d)
        self.vel = np.array(zero3d)
        self.vel_old = np.array(zero3d)
        self.force = np.array(zero3d)
        self.acc = np.array(zero3d)     # acceleration
        
        
        return
        
class sim_config:
    def __init__(self):

        return
class md_system:
    def __init__(self):
        n_freedom = 1
    
class Integrator():
    def __init__(self):
    
        return
        
    def velocityVerlet(self, site, dt):
        """
        Deceptively simple (read about Velocity Verlet on wikipedia)
        Verlet step for one variable
        half step algorithm is used in this version:
        v_{n+0.5} = v_n + 0.5 a_n \Delta t
        x_{n+1} = x_n + v_{n+0.5} \Delta t
        calculate a_{n+1}
        v_{n+1} = v_{n+0.5} + 0.5 a_{n+1} \Delta t
        
        OR eliminate the half step.
        x_{n+1} = x_n + v_n \Delta t + 1/2 a_n \Delta t^2
        and calculate a_{n+1}
        v_{n+1} = v_n + 1/2(a_n + a_{n+1}) \Delta t
        """
        # 3d
        site.vel += 0.5 * site.acc * dt
        site.pos += site.vel * dt
        calc_acc(site from new pos)
        site.vel += 0.5 * site.acc * dt
        
        
        
    def Leapfrog(self, site):
        """
        see http://en.wikipedia.org/wiki/Verlet_integration
        x_{n+1} = 2 x_n - x_{n-1} + a_n \Delta t^2
        when n = 1, we take
        x_1 = x_0 + v_0 \Delta t + 1/2 a_0 \Delta t^2
        with
        v_n = \frac{x_{n+1} - x_n}{\Delta t}
        """
        print "not implement"
        return
        
class BOVDyn():
    def __init__(self):
        """
        some basic structures like structure in c
        """
        self.md = {}
        
        self.md['config'] = {
            "time_step": 0.5,
            
        }
        
        self.md['system'] = {
            "n_site": 2,
            "cluster": []
        }
        site_data = {
                    "pos": [],
                    "vel": [],
                    "force": [],
        }
        
        #
        
        return
        
    def get_md_config(self, keyword = ""):
        """ get the config info """
        cfg = self.md['config']
        if keyword in cfg:
            value = cfg[keyword]
        else:
            print "no such keyword: ", keyword
            exit(1)
        return
    
    
    def update_step_nve(self):
        """
        predict the next step in nve...
        """
        dt = self.get_md_config(keyword = "time_step")
        n_site = self.get_md_system(keyword = "n_site")
        cluster = self.get_md_system(keyword = "cluster")
        # velocity Verlet
        for site in cluster:
            site.vel += 0.5 * site.acc * dt
            site.pos += site.vel * dt
            
        # calculate forces at new pos
        ffmodel.compute_forces(cluster)
        
        # update velocity
        for site in cluster:
            site.vel += 0.5 * site.acc * dt
            
        return
    
    def update_step_ld(self):
        """
        use langevin dynamic method to update md
        see: ..
        """
        
        return
    