#! /usr/bin/env python

# md for rigid body and particles.
#

#! /usr/bin/env python

import os
import math
import random 
import numpy as np

# this is a very very very simple md eigne 
# first @SDU by dulikai
# then  @qibebt dulikai
#
#
# aim at propagation of the nuclear position
#


# struct body {
	# mat_t rotmat;
	# vec_t pos;
	# vec_t vel;
	# vec_t vel_old;
	# vec_t angmom;
	# vec_t angmom_old;
	# vec_t force;
	# vec_t torque;
	# vec_t inertia;
	# vec_t inertia_inv;
	# double mass;
# };

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
        self.angacc = zero3d
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
        zero3d = np.array(3)
        self.pos = np.array(zero3d)
        
        self.md = {}
        
        self.md['config'] = {
            "time_step": 0.5,
            "temperature": 298.15,
            "enable_pbc": "no"
        }
        
        self.md['system'] = {
            "n_sites": 2,
            "sites": [],
            "box": np.zeros(3)
        }
        
        self.md['constant'] = {
            "PI": 3.14159265358979323846,
            "EPSILON": 1.0e-8
        }
        
        self.md['phys'] = {
            "BOLTZMANN": 3.166811429e-6
        # /* Boltzmann constant in [Hartree / K] */
        }
        site_data = {
                    "pos": [],
                    "vel": [],
                    "force": [],
        }
        
        #
        
        return
        
    def line2vars(self, line):
        rec = line.split()
        vel = np.array([0.0, 0.0, 0.0])
        name = rec[0]
        coord = np.array([float(rec[1]), float(rec[2]), float(rec[3])])
        if len(rec) == 7:
            vel = np.array([float(rec[4]), float(rec[5]), float(rec[6])])
        return name, coord, vel
        
        
    def read_coord(self):
        """
        read cart. coordinate in xyz format
        read euler angle in abc format
        """
        sites = []
        # efp xyz
        fp = open("efp.xyz", "r")
        n_sites = int(fp.readline())
        line = fp.readline()
        for i in xrange(n_sites):
            line = fp.readline()
            name, cart, vel = self.line2vars(line)   
            mysite = {'name': name, 'pos': cart, 'vel': vel}
            sites.append(mysite)            
        fp.close()
        # efp abc
        fp = open("efp.abc", "r")
        n_sites = int(fp.readline())
        line = fp.readline()
        for i in xrange(n_sites):
            line = fp.readline()
            name, euler, angvel = self.line2vars(line) 
            sites[i]['angpos'] = euler
            sites[i]['angvel'] = angvel
        self.md['sites'] = sites
        self.md['n_sites'] = n_sites
        self.md['n_freedom'] = n_sites * 6
        fp.close()
        # done
        return
        
    def read_mass(self):
        """
        read in site mass
        """
        sites = self.md['sites']
        # read mass
        fp = open("efp.mass", "r")
        n_sites = int(fp.readline())
        line = fp.readline()
        for i in xrange(n_sites):
            mass = float(fp.readline())
            mysite = sites[i]
            mysite['mass'] = mass           
        fp.close()        
        # read inertia
        fp = open("efp.inertia", "r")
        n_sites = int(fp.readline())
        line = fp.readline()
        for i in xrange(n_sites):
            line = fp.readline()
            rec = line.split()
            inertia = np.array([float(rec[0]), float(rec[1]), float(rec[2])])
            sites[i]['inertia'] = inertia
        fp.close()    
        return
        
        
    def get_md_config(self, keyword=""):
        cfg = self.md['config']
        if keyword in cfg:
            value = cfg[keyword]
        else:
            print("no such md config keyword")
            exit(1)
        return value
        
        
# Initialization 
# initial coordinates: 
# can be obtained from a previous simulation 
# Otherwise: 
# 1. liquid 
# -start out in lattice fcc configuration with appropriate density and random orientation of dipole moments 
# -equilibrate to melt lattice 
# 2. proteins 
# -start with X-ray structure 
# -minimize energy to eliminate any large forces due to repulsive van der Waals contacts or poor geometries 
# -place constraints on all atoms (tethering) to prevent large shifts from starting constraints due to very large forces 
# initial velocities: 
# Gaussian distribution of velocities  ,  , and  appropriate to a given temperature 
# (i.e. Maxwell-Boltzmann distribution) 
# Procedure to generate a Maxwell-Boltzmann distribution of velocities: 
# 1. use a random number generator to generate a normal distribution with zero mean and unit variance (i.e. a set of 3N numbers corresponding to such a normal distribution) 
# 2. multiply by  to get proper distribution of velocities 
# 3. correct so there is no overall momentum 
# There is no unique set of initial velocities - this depends on the random number seed 
# An ensemble of such initial velocity sets satisfies  
        
    def velocitize(self):
        """
        init vel. et al.
        """
        random.seed()
        t = self.get_md_config(keyword="temperature")
        BOLTZMANN = self.md['phys']['BOLTZMANN'] 
        # Boltzmann distribution:
        # every freedom has 1/2 k_B T kinetic energy
        ke = 0.5 * BOLTZMANN * t * 1.0
        sites = self.md['sites']
        for mysite in sites:
            mass = mysite['mass']       
            inertia = mysite['inertia']  
            # 1/2 m v^2 = 1/2 k_B T = k_e
            vel = np.sqrt(2.0 * ke / mass)
            angvel = np.sqrt(2.0 * ke / inertia)
            mysite['vel'] = vel * np.random.normal(0.0, 1.0, 3)
            mysite['angvel'] = angvel * np.random.normal(0.0, 1.0, 3)
        return
        
    def get_kinetic_energy(self):
        """
        T = 0.5 \sum {m v^2}
        """
        n_sites = self.md['n_sites']
        sites = self.md['sites']
        ke = 0.0
        for mysite in sites:
            vel = mysite['vel']
            angvel = mysite['angvel']
            mass = mysite['mass']
            inertia = mysite['inertia']
            vk1 = mass * vel * vel
            vk2 = inertia * angvel * angvel
            ke += np.linalg.norm(vk1)
            ke += np.linalg.norm(vk2)
        return 0.5 * ke
               
    def get_temperature(self):
        """
        return temperature
        """
        BOLTZMANN = self.md['phys']['BOLTZMANN'] 
        n_freedom = self.md['n_freedom']
        ke = self.get_kinetic_energy()
        t = 2.0 * ke / BOLTZMANN / n_freedom
        return t

    def get_dist_mat(self):
        """
        calculate the dist matrix.
        """
        sites = self.md['site']
        n_site = self.md['n_sites']
        self.dist_mat = [[0.0 for j in xrange(n_site)] for i in xrange(n_site)]
        for i in xrange(n_site):
            for j in xrange(n_site):
                coordi = site[i]['coord']
                coordj = site[j]['coord']
                self.dist_mat[i][j] = np.linalg.norm(coordi-coordj)
        
        return
    def get_potential_energy(self):
        """
        Write a function that computes the Lennard-Jones potential
        V = 4 \epsilon [(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6]
        """
        lj.get_distances()
        sigma = 0.3345
        epsilon = 0.0661
        dist_mat = self.dist_mat
        n_site = self.cluster['n_site']
        v = 0.0
        for i in xrange(n_site):
            for j in xrange(i+1, n_site):
                r = dist_mat[i][j]
                s = sigma / r
                s6 = s**6; s12 = s6*s6;
                v += 4.0 * epsilon * (s12 - s6)
        return v
        
    def get_energy_nve(self):
        """
        return total energy
        """
        pot = self.get_potential_energy()
        ke = self.get_kinetic_energy()        
        return pot + ke
        
    def boxize(self, mysite):
        """
        enable pbc and wrap the 
        """
        pbc = self.get_md_config(keyword="enable_pbc")
        box = self.md['system']['box']
        pos = mysite['pos']
        if pbc == "yes":
            mypos = pos
        else:
            mypos = pos + box.astype(int) * box
        return mypos
        

    def get_system_com(self):
        mass = 0.0
        com = np.zeros(3)
        sites = self.md['sites']
        for mysite in sites:
            pos = self.boxize(mysite)
            com += mysite['mass'] * pos
            mass += mysite['mass']
        com /= mass
        return com
            
    def get_system_com_velocity(self):
        """
        velocity of the system, trans.
        """
        tmass = 0.0
        cov = np.zeros(3)
        sites = self.md['sites']
        for mysite in sites:
            vel = mysite['vel']
            mass = mysite['mass']
            cov = vel * mass
            tmass += mass
        cov /= tmass
        return cov
    
    def get_system_angvel(self):
        """
        angular velocity of the system
        """
        angvel = np.zeros(3)
        sites = self.md['sites']
        for mysite in sites:
            angvel += mysite['angvel']
        return angvel

    def get_system_angmom(self):
        """
        get system angular moment rot.
        """
        com = self.get_system_com()
        cov = self.get_system_velocity()
        tangmom = np.zeros(3)
        sites = self.md['sites']
        for mysite in sites:
            pos = self.boxize(mysite)
            mass = mysite['mass']
            dr = pos - com
            dv = vel - cov
            angmom += np.cross(dr, dv) * mass
        return angmom

    def get_inertia_tensor(self, vr, mass):
        """
        the moment of inertia for one mass piece
        """
        x, y, z = vr
        xx = mass * (y * y + z * z)
        yy = mass * (x * x + z * z)
        zz = mass * (x * x + y * y)
        xy = - mass * x * y
        xz = - mass * x * z
        yz = - mass * y * z
        yx = xy
        xz = zx
        yz = zy
        inertia = [
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ]
        return inertia
        
    def get_system_inertia_tensor(self):
        """
        system inertia ..
        """
        inertia = np.zeros((3,3))
        com = self.get_system_com()
        sites = self.md['sites']
        for mysite in sites:
            pos = self.boxize(mysite)
            mass = mysite['mass']
            dr = pos - com
            inertia += self.get_inertia_tensor(dr, mass)
        return
    

    def remove_system_drift(self):
        """
        remove the drift error
        """
        com = self.get_system_com()
        cov = self.get_system_com_velocity()
        angvel = self.get_system_angvel()
        EPSILON = 1.0e-8
        inertia = self.get_system_inertia_tensor()
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
        
    def update_step_nve(self):
        """
        predict the next step in nve...
        for rigid body, 
        propagate the euler angle & euler velocity (angular vel.)
        """
        dt = self.get_md_config(keyword = "time_step")
        cluster = self.get_md_system(keyword = "cluster")
        # velocity Verlet
        for site in cluster:
            site.vel += 0.5 * site.acc * dt
            site.angvel += 0.5 * site.angacc * dt
            site.pos += site.vel * dt
            site.angpos += site.angvel * dt
        # calculate forces at new pos
        ffmodel.compute_forces(cluster, type="simple")        
        # update velocity
        for site in cluster:
            site.vel += 0.5 * site.acc * dt
            site.angvel += 0.5 * site.angacc * dt
            
        return
        
    
    def velocityVerlet(self):
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
        
    def update_step_nvt(self):
        """
        NVT with Nose-Hoover thermostat:
        William G. Hoover
        Canonical dynamics: Equilibrium phase-space distributions
        Phys. Rev. A 31, 1695 (1985)
        """
        dt = self.get_md_config(keyword = "time_step")
        target = self.get_md_config(keyword = "temperature")
        tau = self.get_md_config(keyword = "thermostat_tau")
        chi = self.get_md_config(keyword = "nvt_chi")
        chi_dt = self.get_md_config(keyword = "nvt_chi_dt")
        t0 = self.get_temperature()
        sites = self.md['sites']
        # velocity Verlet
        for site in cluster:
            site.vel += 0.5 * dt * (site.acc - site.vel * chi)
            site.angvel += 0.5 * site.angacc * dt
            site.pos += site.vel * dt
            site.angpos += site.angvel * dt
            
		body->vel.x += 0.5 * dt * (body->force.x / body->mass -
					body->vel.x * data->chi);
		body->vel.y += 0.5 * dt * (body->force.y / body->mass -
					body->vel.y * data->chi);
		body->vel.z += 0.5 * dt * (body->force.z / body->mass -
					body->vel.z * data->chi);

		body->angmom.x += 0.5 * dt * (body->torque.x -
					body->angmom.x * data->chi);
		body->angmom.y += 0.5 * dt * (body->torque.y -
					body->angmom.y * data->chi);
		body->angmom.z += 0.5 * dt * (body->torque.z -
					body->angmom.z * data->chi);

		body->pos.x += body->vel.x * dt;
		body->pos.y += body->vel.y * dt;
		body->pos.z += body->vel.z * dt;

		rotate_body(body, dt);
            
        # calculate forces at new pos
        ffmodel.compute_forces(cluster, type="simple")        
        # update velocity
        for site in cluster:
            site.vel += 0.5 * site.acc * dt
            site.angvel += 0.5 * site.angacc * dt
        

/*
 * NVT with Nose-Hoover thermostat:
 *
 * William G. Hoover
 *
 * Canonical dynamics: Equilibrium phase-space distributions
 *
 * Phys. Rev. A 31, 1695 (1985)
 */
static void update_step_nvt(struct md *md)
{
	struct nvt_data *data = (struct nvt_data *)md->data;

	double dt = cfg_get_double(md->state->cfg, "time_step");
	double target = cfg_get_double(md->state->cfg, "temperature");
	double tau = cfg_get_double(md->state->cfg, "thermostat_tau");

	double t0 = get_temperature(md);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * dt * (body->force.x / body->mass -
					body->vel.x * data->chi);
		body->vel.y += 0.5 * dt * (body->force.y / body->mass -
					body->vel.y * data->chi);
		body->vel.z += 0.5 * dt * (body->force.z / body->mass -
					body->vel.z * data->chi);

		body->angmom.x += 0.5 * dt * (body->torque.x -
					body->angmom.x * data->chi);
		body->angmom.y += 0.5 * dt * (body->torque.y -
					body->angmom.y * data->chi);
		body->angmom.z += 0.5 * dt * (body->torque.z -
					body->angmom.z * data->chi);

		body->pos.x += body->vel.x * dt;
		body->pos.y += body->vel.y * dt;
		body->pos.z += body->vel.z * dt;

		rotate_body(body, dt);
	}

	data->chi += 0.5 * dt * (t0 / target - 1.0) / tau / tau;
	data->chi_dt += 0.5 * dt * data->chi;

	compute_forces(md);

	double chi_init = data->chi;
	vec_t angmom_init[md->n_bodies], vel_init[md->n_bodies];

	for (size_t i = 0; i < md->n_bodies; i++) {
		angmom_init[i] = md->bodies[i].angmom;
		vel_init[i] = md->bodies[i].vel;
	}

	for (size_t iter = 1; iter <= MAX_ITER; iter++) {
		double chi_prev = data->chi;
		double ratio = get_temperature(md) / target;

		data->chi = chi_init + 0.5 * dt * (ratio - 1.0) / tau / tau;

		for (size_t i = 0; i < md->n_bodies; i++) {
			struct body *body = md->bodies + i;

			body->vel.x = vel_init[i].x + 0.5 * dt *
				(body->force.x / body->mass - vel_init[i].x * data->chi);
			body->vel.y = vel_init[i].y + 0.5 * dt *
				(body->force.y / body->mass - vel_init[i].y * data->chi);
			body->vel.z = vel_init[i].z + 0.5 * dt *
				(body->force.z / body->mass - vel_init[i].z * data->chi);

			body->angmom.x = angmom_init[i].x + 0.5 * dt *
				(body->torque.x - angmom_init[i].x * data->chi);
			body->angmom.y = angmom_init[i].y + 0.5 * dt *
				(body->torque.y - angmom_init[i].y * data->chi);
			body->angmom.z = angmom_init[i].z + 0.5 * dt *
				(body->torque.z - angmom_init[i].z * data->chi);
		}

		if (fabs(data->chi - chi_prev) < EPSILON)
			break;

		if (iter == MAX_ITER)
			msg("WARNING: NVT UPDATE DID NOT CONVERGE\n\n");
	}

	data->chi_dt += 0.5 * dt * data->chi;
}


    def print_info(self):
        """
        print dynamics info.
        """
        pot = self.md['potential_energy']
        kin = self.get_kinetic_energy()
        ene = self.get_energy_nve()
        temperature = self.get_temperature()        
        print("%30s %16.10lf\n", "POTENTIAL ENERGY", pot)
        print("%30s %16.10lf\n", "KINETIC ENERGY", kin)
        print("%30s %16.10lf\n", "INVARIANT", ene)
        print("%30s %16.10lf\n", "TEMPERATURE", temperature)
        return
        
    def print_restart(self):
        """
        restart dump
        """
        

static void print_restart(const struct md *md)
{
	msg("    RESTART DATA\n\n");

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		char name[64];
		check_fail(efp_get_frag_name(md->state->efp, i, sizeof(name), name));

		double xyzabc[6] = { body->pos.x * BOHR_RADIUS,
				     body->pos.y * BOHR_RADIUS,
				     body->pos.z * BOHR_RADIUS };

		matrix_to_euler(&body->rotmat, xyzabc + 3, xyzabc + 4, xyzabc + 5);

		double vel[6] = { body->vel.x,
				  body->vel.y,
				  body->vel.z,
				  body->angmom.x * body->inertia_inv.x,
				  body->angmom.y * body->inertia_inv.y,
				  body->angmom.z * body->inertia_inv.z };

		print_fragment(name, xyzabc, vel);
	}

	msg("\n");
}

static struct md *md_create(struct state *state)
{
	struct md *md = xcalloc(1, sizeof(struct md));

	md->state = state;
	md->box = box_from_str(cfg_get_string(state->cfg, "periodic_box"));

	switch (cfg_get_enum(state->cfg, "ensemble")) {
		case ENSEMBLE_TYPE_NVE:
			md->get_invariant = get_invariant_nve;
			md->update_step = update_step_nve;
			break;
		case ENSEMBLE_TYPE_NVT:
			md->get_invariant = get_invariant_nvt;
			md->update_step = update_step_nvt;
			md->data = xcalloc(1, sizeof(struct nvt_data));
			break;
		case ENSEMBLE_TYPE_NPT:
			md->get_invariant = get_invariant_npt;
			md->update_step = update_step_npt;
			md->data = xcalloc(1, sizeof(struct npt_data));
			break;
		default:
			assert(0);
	}

	md->n_bodies = state->sys->n_frags;
	md->bodies = xcalloc(md->n_bodies, sizeof(struct body));

	double coord[6 * md->n_bodies];
	check_fail(efp_get_coordinates(state->efp, coord));

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->pos.x = coord[6 * i + 0];
		body->pos.y = coord[6 * i + 1];
		body->pos.z = coord[6 * i + 2];

		double a = coord[6 * i + 3];
		double b = coord[6 * i + 4];
		double c = coord[6 * i + 5];

		euler_to_matrix(a, b, c, &body->rotmat);

		body->vel.x = md->state->sys->frags[i].vel[0];
		body->vel.y = md->state->sys->frags[i].vel[1];
		body->vel.z = md->state->sys->frags[i].vel[2];

		set_body_mass_and_inertia(state->efp, i, body);

		body->angmom.x = md->state->sys->frags[i].vel[3] * body->inertia.x;
		body->angmom.y = md->state->sys->frags[i].vel[4] * body->inertia.y;
		body->angmom.z = md->state->sys->frags[i].vel[5] * body->inertia.z;

		md->n_freedom += 3;

		if (body->inertia.x > EPSILON)
			md->n_freedom++;
		if (body->inertia.y > EPSILON)
			md->n_freedom++;
		if (body->inertia.z > EPSILON)
			md->n_freedom++;
	}

	return (md);
}

static void print_status(const struct md *md)
{
	print_geometry(md->state->efp);
	print_restart(md);
	print_info(md);

	fflush(stdout);
}

static void md_shutdown(struct md *md)
{
	free(md->bodies);
	free(md->data);
	free(md);
}

void print_geometry(struct efp *efp)
{
	size_t n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	msg("    GEOMETRY (ANGSTROMS)\n\n");

	for (size_t i = 0; i < n_frags; i++) {
		size_t n_atoms;
		check_fail(efp_get_frag_atom_count(efp, i, &n_atoms));

		struct efp_atom atoms[n_atoms];
		check_fail(efp_get_frag_atoms(efp, i, n_atoms, atoms));

		for (size_t a = 0; a < n_atoms; a++) {
			double x = atoms[a].x * BOHR_RADIUS;
			double y = atoms[a].y * BOHR_RADIUS;
			double z = atoms[a].z * BOHR_RADIUS;

			msg("%-16s %12.6lf %12.6lf %12.6lf\n", atoms[a].label, x, y, z);
		}
	}

	size_t n_charges;
	check_fail(efp_get_point_charge_count(efp, &n_charges));

	if (n_charges > 0) {
		double xyz[3 * n_charges];
		check_fail(efp_get_point_charge_coordinates(efp, xyz));

		for (size_t i = 0; i < n_charges; i++) {
			char label[32];
			double x = xyz[3 * i + 0] * BOHR_RADIUS;
			double y = xyz[3 * i + 1] * BOHR_RADIUS;
			double z = xyz[3 * i + 2] * BOHR_RADIUS;

			snprintf(label, sizeof(label), "Q%04zu", i + 1);
			msg("%-16s %12.6lf %12.6lf %12.6lf\n", label, x, y, z);
		}
	}

	msg("\n\n");
}

void sim_md(struct state *state)
{
	msg("MOLECULAR DYNAMICS JOB\n\n\n");

	struct md *md = md_create(state);

	if (cfg_get_bool(state->cfg, "velocitize"))
		velocitize(md);

	remove_system_drift(md);
	compute_forces(md);

	msg("    INITIAL STATE\n\n");
	print_status(md);

	for (int i = 1; i <= cfg_get_int(state->cfg, "max_steps"); i++) {
		md->update_step(md);

		if (i % cfg_get_int(state->cfg, "print_step") == 0) {
			msg("    STATE AFTER %d STEPS\n\n", i);
			print_status(md);
		}
	}

	md_shutdown(md);

	msg("MOLECULAR DYNAMICS JOB COMPLETED SUCCESSFULLY\n");
}
        
    # def get_md_config(self, keyword = ""):
        # """ get the config info """
        # cfg = self.md['config']
        # if keyword in cfg:
            # value = cfg[keyword]
        # else:
            # print "no such keyword: ", keyword
            # exit(1)
        # return

       

        
        
    # state(md)    
    # velocitize(md);
    # remove_system_drift(md)
	# compute_forces(md);
    
    def sim_md(self):
        """
        MOLECULE DYNAMICS JOBS
        """
        # build working state
        self.setup()
        # create md vars
        self.create()
        # initial velocity
        self.velocitize()
        # remove system drift
        self.remove_drift()
        # then, the force
        ffmodel.compute_forces(self.md, type="simple")
        #
        # start ...
        print "INITIAL STATE"
        self.print_status()
        
        n_step
        for i in xrange(n_step):
            self.update_step()
            if i % n_print_step == 0:
                print "    STATE AFTER %d STEPS\n\n"
                self.print_status()
                
        self.md_shutdown()
        print "MOLECULAR DYNAMICS JOB COMPLETED SUCCESSFULLY\n"
        return 
        
        
# void sim_md(struct state *state)
# {
	# msg("MOLECULAR DYNAMICS JOB\n\n\n");

	# struct md *md = md_create(state);

	# if (cfg_get_bool(state->cfg, "velocitize"))
		# velocitize(md);

	# remove_system_drift(md);
	# compute_forces(md);

	# msg("    INITIAL STATE\n\n");
	# print_status(md);

	# for (int i = 1; i <= cfg_get_int(state->cfg, "max_steps"); i++) {
		# md->update_step(md);

		# if (i % cfg_get_int(state->cfg, "print_step") == 0) {
			# msg("    STATE AFTER %d STEPS\n\n", i);
			# print_status(md);
		# }
	# }

	# md_shutdown(md);

	# msg("MOLECULAR DYNAMICS JOB COMPLETED SUCCESSFULLY\n");
# }
        
    def update_step_ld(self):
        """
        use langevin dynamic method to update md
        see: ..
        """
        
        return
    
    
if __name__ == "__main__":
    dyn = BOVDyn()
    dyn.read_coord()
    dyn.read_mass()
    dyn.velocitize()
    t = dyn.get_temperature()
    
    
    print dyn.md
    print t
    
    
    
    
    