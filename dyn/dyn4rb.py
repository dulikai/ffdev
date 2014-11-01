#! /usr/bin/env python

# md for rigid body and particles.
#

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
        zero3d = np.array([0.0, 0.0, 0.0])
        self.pos = np.array(zero3d)
        self.pos_old = np.array(zero3d)
        self.vel = np.array(zero3d)
        self.vel_old = np.array(zero3d)
        self.force = np.array(zero3d)
        self.acc = np.array(zero3d)     # acceleration
        
        self.angmom = np.array(zero3d)
        self.angvel = np.array(zero3d)
        self.ang = zero3d
        
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
static void update_step_nve(struct md *md)
{
	double dt = cfg_get_double(md->state->cfg, "time_step");

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angmom.x += 0.5 * body->torque.x * dt;
		body->angmom.y += 0.5 * body->torque.y * dt;
		body->angmom.z += 0.5 * body->torque.z * dt;

		body->pos.x += body->vel.x * dt;
		body->pos.y += body->vel.y * dt;
		body->pos.z += body->vel.z * dt;

		rotate_body(body, dt);
	}

	compute_forces(md);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angmom.x += 0.5 * body->torque.x * dt;
		body->angmom.y += 0.5 * body->torque.y * dt;
		body->angmom.z += 0.5 * body->torque.z * dt;
	}
}

    def rotate_body(self):
        
    
    def update_step_nve_simple(self):
        """
        predict the next step in nve...
        for rigid body, 
        propagate the euler angle & euler velocity (angular vel.)
        """
        dt = self.get_md_config(keyword = "time_step")
        n_site = self.get_md_system(keyword = "n_site")
        cluster = self.get_md_system(keyword = "cluster")
        # velocity Verlet
        for site in cluster:
            site.vel += 0.5 * site.acc * dt
            site.angmom += 0.5 * site.torque * dt
            site.pos += site.vel * dt
            self.rotate_body()
            
        # calculate forces at new pos
        ffmodel.compute_forces(cluster)
        
        # update velocity
        for site in cluster:
            site.vel += 0.5 * site.acc * dt
            
        return
    
    
    def remove_system_drift2(self):
        """
        remove drift error angular momentum
        """
        com = self.get_system_com()
        cov = self.get_system_com_velocity()
        angmom = self.get_system_angmom()
        inertia_inv = np.zeros((3,3))
        EPSILON = 1.0e-8
        inertia = self.get_system_inertia_tensor()
        det = np.linalg.det(inertia)
        if det < EPSILON:
            print "Inertia det too small... <continue..>"
        inertia_inv = np.linalg.inv(inertia)
        angvel = np.dot(inertia_inv, angmom)
        sites = self.md['sites']
        for mysite in sites:
            pos = self.boxize(mysite)
            mass = mysite['mass']
            dr = pos - com
            dv = np.cross(angvel, dr)
            mysite['vel'] -= cov + dv
        cov2 = self.get_system_com_velocity()
        angvel2 = self.get_system_angvel()
        if np.linalg.norm(cov2) > EPSILON || np.linalg.norm(angvel2) > EPSILON:
            exit(1)
        return
    def velocityVerlet(self):
        """
        useless
        """
        dt = 0.1
        for i in xrange(1000):
            energy = self.interaction()
            self.compute_nd_acc()
            print "Step %d E=%20.12f" % (i, energy)
            acc_old = []
            self.update_step_nve()
            for i_site in xrange(n_site):
                vel = site[i_site]['vel']
                acc = site[i_site]['acc']
                coord = site[i_site]['coord']
                print coord    
                coord += vel * dt + acc * dt * dt * 0.5
                # print site[i_site]['coord']            
                acc_old.append(acc.copy())
                
            self.compute_nd_acc(sigma, epsilon)
            for i_site in xrange(n_site):
                vel = site[i_site]['vel']
                acc = site[i_site]['acc']
                vel += (acc_old[i_site] + acc) * dt * 0.5 
                # print acc_old-acc
        return  



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
        