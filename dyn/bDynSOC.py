#! /usr/bin/env python
#
#
# 2014.11
# this is a very very very simple md eigne for SOC style hopping method
# first @SDU by dulikai 2012
# then  @qibebt dulikai 2014
#
# aim at propagation of the Newton-style equation with surface hopping 
# of singlet triplet etc.
# md for particles.
#
# http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/Newtonhtml.html
#
#

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
# Gaussian distribution of velocities, appropriate to a given temperature 
# Maxwell-Boltzmann distribution):
# Procedure to generate a Maxwell-Boltzmann distribution of velocities: 
# 1. use a random number generator to generate a normal distribution 
# with zero mean and unit variance 
# (i.e. a set of 3N numbers corresponding to such a normal distribution) 
# 2. multiply by  to get proper distribution of velocities 
# 3. correct, so there is no overall momentum 
# There is no unique set of initial velocities - this depends on the random number seed 
# An ensemble of such initial velocity sets satisfies  

#
#
import os
import math
import random 
import numpy as np
import copy

from bStatus import bStatus
from bConfig import bConfig
from bConst import bConst
from ffTest import ffTest

from quantum.gaussian import Gaussian
from quantum import tools

class bDynSOC():
    """
    this model implement a very-very simple md runner
    """
    # internal data set.
    data = {
        "nvt_chi": 0.0,
        "nvt_chi_dt": 0.0,
        "n_max_nvt": 10,
    }
    
    def __init__(self, config):
        """
        some basic structures like structure in c
        """
        self.config = config    # md control
        # self.model = model      # dynamic variable & ff parameters
        return

    def prob_isc(self, Hso, gdiff, vel):
        """
        Landau Zener model 
        """
        t = np.dot(gdiff, vel)
        xi = 8.0 * Hso * Hso / t
        prob = 1.0 - np.exp(-np.pi / 4.0 * xi)
        return prob
        
    def hopping(self, pes, grad, Hso):
        """
        calc. hopping prob.
        """
        i_state = 1
        for i in xrange(n_state):
            pass
        
        # Reject all hops when the energy gap is large than hop_e if necessary
        if (abs(pes_all(old_index_state) - pes_all(new_index_state)) >= hop_energy):
            new_index_state = old_index_state
            index_state = old_index_state

        return

    def setup(self):
        """
        initial based on configure
        """
        os.chdir("test")
        print "I am here, ", os.getcwd()
        self.status = bStatus()
        np.random.seed()
        return

    def dump(self):
        sites = self.status.give(keyword="sites")
        n_site = self.status.give(keyword = "n_site")
        parm = self.config.give(keyword="quantum")
        parm['n_atom'] = n_site
        mol = {}
        mol['natom'] = n_site
        atoms = []
        for mysite in sites:
            name = mysite.name
            coord = list(mysite.pos)
            atoms.append({'name': name, 'coord': coord})
        mol['atoms'] = atoms
        obj = {'parm': parm, 'mol': mol}
        tools.dump_data("interface.json", obj)
        print "DUMP INTERFACE FILE"
        return
            
    def aux_action(self):
        """
        based on force, update acc. et al.
        """
        # other case, need to map model <--> status
        sites = self.status.give(keyword="sites")
        # self.status.potential_energy = 1.0  # set as default...
        for mysite in sites:
            # mysite.force = np.random.normal(0.0, 1.0, 3)
            mysite.acc = mysite.force / mysite.mass
        return
        
    def action(self):
        """ the action by potential energy """
        # dump interface file...
        self.dump()
        # call the energy & gradient worker
        ff = ffTest(self.status)
        ff.eandg()
        return

    def velocitize(self):
        """
        init vel. et al.
        """
        np.random.seed()
        Boltzmann = bConst.k_B
        temperature = self.config.give(keyword="temperature")
        n_freedom = self.status.give(keyword = "n_freedom")
        n_site = self.status.give(keyword = "n_site")
        sites = self.status.give(keyword="sites")
        # Boltzmann distribution:
        # every freedom has 1/2 k_B T = 1/2 m <v^2> kinetic energy
        ke = 0.5 * Boltzmann * temperature #* n_freedom / (6.0 * n_site)
        for mysite in sites:
            mass = mysite.mass       
            inertia = mysite.inertia  
            # 1/2 m v^2 = 1/2 k_B T = k_e
            vel = np.sqrt(2.0 * ke / mass)
            angvel = np.sqrt(2.0 * ke / inertia)
            # angmom = np.sqrt(2.0 * ke * inertia)
            mysite.vel = vel * np.random.normal(0.0, 1.0, 3)
            mysite.angvel = angvel * np.random.normal(0.0, 1.0, 3)
            mysite.angmom = inertia * mysite.angvel
        return
        
    def get_kinetic_energy(self):
        """
        T = 0.5 \sum {m v^2}
        """
        sites = self.status.give(keyword="sites")
        kinetic = 0.0
        for mysite in sites:
            vel = mysite.vel
            angvel = mysite.angvel
            mass = mysite.mass
            inertia = mysite.inertia
            vk1 = mass * vel * vel
            vk2 = inertia * angvel * angvel
            kinetic += np.linalg.norm(vk1)
            kinetic += np.linalg.norm(vk2)
        #print kinetic
        return 0.5 * kinetic
               
    def get_temperature(self):
        """
        return temperature
        """
        Boltzmann = bConst.k_B
        n_freedom = self.status.give(keyword = "n_freedom")
        kinetic = self.get_kinetic_energy()
        t = 2.0 * kinetic / Boltzmann / n_freedom
        return t

    def get_dist_mat(self):
        """
        calculate the distance matrix.
        """
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        dist_mat = [[0.0 for j in xrange(n_site)] for i in xrange(n_site)]
        for i in xrange(n_site):
            for j in xrange(n_site):
                ri = sites[i].pos
                rj = sites[j].pos
                dist_mat[i][j] = np.linalg.norm(ri-rj)
                # print ri, rj
        return dist_mat
        
    def get_potential_energy(self):
        """
        Write a function that computes the Lennard-Jones potential
        V = 4 \epsilon [(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6]
        """
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        sigma = 0.3345; epsilon = 0.0661
        dist = self.get_dist_mat()
        # print dist
        v = 0.0
        for i in xrange(n_site):
            for j in xrange(i+1, n_site):
                r = dist[i][j]
                s = sigma / r
                s6 = s**6; s12 = s6*s6;
                v += 4.0 * epsilon * (s12 - s6)
        return v
        
    def get_energy_nve(self):
        """
        return total energy
        """
        potential_energy = self.get_potential_energy()
        kinetic_energy = self.get_kinetic_energy()        
        return potential_energy + kinetic_energy
        
    def boxize(self, mysite):
        """
        enable pbc and wrap it
        """
        mypos = mysite.pos
        pbc = self.config.give(keyword="enable_pbc")
        if pbc == "yes":
            box = self.status.give(keyword="box")
            mypos += (pos/box).astype(int) * box
        return mypos
        
    def get_system_com(self):
        """ center of mass """
        sites = self.status.give(keyword="sites")
        com = np.zeros(3)
        mass = 0.0
        for mysite in sites:
            pos = self.boxize(mysite)
            m = mysite.mass
            com += m * pos
            mass += m
        com /= mass
        return com
            
    def get_system_cov(self):
        """
        com velocity of the system, trans.
        """
        sites = self.status.give(keyword="sites")
        cov = np.zeros(3)
        mass = 0.0
        for mysite in sites:
            v = mysite.vel
            m = mysite.mass
            cov += v * m
            mass += m
        cov /= m
        return cov
    
    def get_inertia_tensor_contrib(self, vr, mass):
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
        zx = xz
        zy = yz
        inertia = np.array([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ])
        return inertia
        
    def get_system_inertia_tensor(self):
        """
        system inertia ..
        """
        com = self.get_system_com()
        sites = self.status.give(keyword="sites")
        inertia = np.zeros((3,3))
        for mysite in sites:
            pos = self.boxize(mysite)
            dr = pos - com
            # inertia..
            inertia += self.get_inertia_tensor_contrib(dr, mysite.mass)
        return inertia
    
    def get_system_angmom(self):
        """
        get system angular moment rot.
        """
        com = self.get_system_com()
        cov = self.get_system_cov()
        sites = self.status.give(keyword="sites")
        angmom = np.zeros(3)
        for mysite in sites:
            pos = self.boxize(mysite)
            vel = mysite.vel
            dr = pos - com
            dv = vel - cov
            angmom += np.cross(dr, dv) * mysite.mass
        return angmom

    def get_system_angvel(self):
        """
        angular velocity of the system
        """
        EPSILON = bConst.EPSILON
        com = self.get_system_com()
        cov = self.get_system_cov()
        angmom = self.get_system_angmom()        
        inertia = self.get_system_inertia_tensor()
        
        inertia_inv = np.zeros((3,3))
        if np.linalg.det(inertia) < EPSILON:
            inertia_inv[0][0] = 0.0 if inertia[0][0] < EPSILON else 1.0 / inertia[0][0]
            inertia_inv[1][1] = 0.0 if inertia[1][1] < EPSILON else 1.0 / inertia[1][1]
            inertia_inv[2][2] = 0.0 if inertia[2][2] < EPSILON else 1.0 / inertia[2][2]
        else:
            inertia_inv = np.linalg.inv(inertia)
        angvel = np.dot(inertia_inv, angmom)
        return angvel
                
    def remove_system_drift(self):
        """
        remove the drift error
        """
        com = self.get_system_com()
        cov = self.get_system_cov()
        angvel = self.get_system_angvel()
        angmom = self.get_system_angmom()

        EPSILON = bConst.EPSILON
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        for mysite in sites:
            pos = self.boxize(mysite)    
            dv = np.cross(angvel, pos - com)  
            mysite.vel -= cov / n_site + dv
        cov2 = self.get_system_cov()
        angmom2 = self.get_system_angmom()
        
        if np.linalg.norm(cov2) < EPSILON and np.linalg.norm(angmom2) < EPSILON:
            print("@note: system drift has been removed successfully")
        else:
            print("@warning: system drift cannot be removed, ...")

        return

    def update_step_nve(self):
        """
        ##
        Deceptively simple (read about Velocity Verlet on wikipedia)
        Verlet update step
        >>half step algorithm [used in this implementation]:
        v_{n+0.5} = v_n + 0.5 a_n \Delta t
        x_{n+1} = x_n + v_{n+0.5} \Delta t
        calculate a_{n+1}
        v_{n+1} = v_{n+0.5} + 0.5 a_{n+1} \Delta t
        
        >>OR eliminate the half step.
        x_{n+1} = x_n + v_n \Delta t + 1/2 a_n \Delta t^2
        and calculate a_{n+1}
        v_{n+1} = v_n + 1/2(a_n + a_{n+1}) \Delta t
        ------
        predict the next step in NVE ensemble
        for rigid body, 
        propagate the euler angle & angular velocity (angular vel.)
        """
        dt = self.config.give(keyword="time_step")
        node = self.status.give(keyword="sites")    
        # velocity Verlet
        for site in node:
            site.vel += 0.5 * site.acc * dt
            site.angvel += 0.5 * site.angacc * dt
            site.pos += site.vel * dt
            site.angpos += site.angvel * dt 
            site.angpos[0] %= (bConst.PI * 2) # alpha [-pi,pi]
            site.angpos[1] %= (bConst.PI * 1) # beta  [0,pi] or [-pi/2,pi/2]
            site.angpos[2] %= (bConst.PI * 2) # gamma [-pi,pi]
            # print "nac", site.force
        # calculate forces at new pos
        self.action()        
        # update velocity
        for site in node:
            site.vel += 0.5 * site.acc * dt
            site.angvel += 0.5 * site.angacc * dt            
        return
        
    def langevin_modify(self):
        """
        update langevin to modify the obtained force & torque
        http://en.wikipedia.org/wiki/Langevin_dynamics
        see: http://en.wikipedia.org/wiki/Stokes%27_law
        http://scitation.aip.org/content/aip/journal/jcp/130/23/10.1063/1.3149788
        THE JOURNAL OF CHEMICAL PHYSICS 130, 234101 2009
        Langevin dynamics for rigid bodies of arbitrary shape
        J. Chem. Phys. 128, 234107 (2008); http://dx.doi.org/10.1063/1.2936991
        or
        http://hyperphysics.phy-astr.gsu.edu/hbase/airfri.html
        http://lammps.sandia.gov/doc/fix_viscous.html
        If an object moves slowly through a liquid, 
        it will experience a viscous drag 
        which is approximately proportional to the velocity. 
        For a sphere falling through a liquid in a velocity regime 
        where viscosity is the dominant effect and 
        turbulence can be neglected.
        the relationship between viscosity (\eta) & friction (\gamma):
        Stokes' law
        \gamma = 6 \pi a \eta
        where a is sphere radius, \eta is the dynamic viscosity of the frictional fluid
        or
        \gamma = k_b T / D
        where D = \frac{k_bT}{6 \pi \eta a}, which is particle diffusion coefficient.
        -----
        \gamma is in the unit of force/velocity, \eta is Pa S
        .....
        Langevin equation
        $$m \ddot{r} = -\Delta V - \gamma m \dot{r} + \sqrt{2 \gamma k_b T m} R(t) $$
        detals: http://en.wikipedia.org/wiki/Langevin_dynamics
        notes: in this version, gamma is constant for every particle
        """
        kb = bConst.k_B
        au2k = bConst.au2k  # a.u. conver to k constant
        gamma = self.config.give(keyword="ld_gamma")
        temperature = self.config.give(keyword="ld_temperature")
        n_site = self.status.give(keyword="n_site")
        sites = self.status.give(keyword="sites")
        for body in sites:
            # trans.
            ran = np.random.normal(3)
            # - \gamma m \dot{r}
            body.force -= body.vel * body.mass * gamma
            # \sqrt{2 \gamma k_b T m} R(t)
            t = 2.0 * kb * (temperature/au2k)
            body.force += np.sqrt(gamma * body.mass * t) * ran
            # update
            body.acc = body.force / body.mass
            # rot.
            ran = np.random.normal(3)
            # \tao = -\Delta V - \gamma I \omega + ...
            body.torque -= body.angvel * body.inertia * gamma
            body.torque += np.sqrt(gamma * body.inertia * t) * ran
            # update
            body.angacc = body.torque / body.inertia
        return
                        
    def update_step_ld(self):
        """
        use langevin dynamics method to update md: NVT 
        units:
        time: ps; 
        temperature: K; 
        nm for length
        ps for time
        amu (g/mol) for mass
        the charge of a proton for charge
        rad for angles
        k_B = 0.0019872041(18) kcal/mol/K 
        or
        kB = 0.0083144621(75 # Boltzmann's constant (kJ/mol/K).
        per mole form often used in statistical mechanics 
        using thermochemical calorie = 4.184 Joule 
        http://en.wikipedia.org/wiki/Boltzmann_constant        
        ----
        """
        sites = self.status.give(keyword="sites")    
        dt = self.config.give(keyword="time_step")
        self.langevin_modify()
        for body in sites:
            body.vel += 0.5 * body.acc * dt
            body.angvel += 0.5 * body.angacc * dt
            body.pos += body.vel * dt
            body.angpos += body.angvel * dt
        # calculate forces at new pos
        self.action()  
        self.langevin_modify()
        # update velocity
        for body in sites:
            body.vel += 0.5 * body.acc * dt
            body.angvel += 0.5 * body.angacc * dt
        return
                
    def update_step_nvt(self):
        """
        NVT with Nose-Hoover thermostat:
        William G. Hoover
        Canonical dynamics: Equilibrium phase-space distributions
        Phys. Rev. A 31, 1695 (1985)
        formula: Berendsen and Nose-Hoover thermostats
        accerlation: \ddot{r} = F_i/m_i  - \chi \dot{r} 
        \tao^2 = \frac{Q}{N k_B T_0} (g/N) : (Nose-Hoover formalism; g=N)
        \dot{\chi} = 1/\tao^2 (T/T_0 - 1)
        site layer page node cluster group 
        """
        EPSILON = bConst.EPSILON
        dt = self.config.give(keyword="time_step")
        target = self.config.give(keyword="temperature")
        tau = self.config.give(keyword="thermostat_tau")
        t0 = self.get_temperature()
        chi = self.data['nvt_chi']
        chi_dt = self.data['nvt_chi_dt']
        n_max = self.data['n_max_nvt']
        node = self.status.give(keyword="sites")    
        # velocity Verlet
        for site in node:
            site.vel += 0.5 * dt * (site.acc - site.vel * chi)
            site.angvel += 0.5 * dt * (site.angacc - site.angvel * chi)
            site.pos += site.vel * dt
            site.angpos += site.angvel * dt
        # propagate \chi
        chi += 0.5 * dt * (t0 / target - 1.0) / tau / tau
        chi_dt += 0.5 * dt * chi
        # calculate forces at new pos
        self.action() 
        # 
        chi_init = chi
        node_init = copy.deepcopy(node)
        # update velocity
        for i in xrange(n_max):
            chi_old = chi
            ratio = self.get_temperature() / target
            chi = chi_init + 0.5 * dt * (ratio - 1.0) / tau / tau;
            for site, site_init in zip(node, node_init):
                site.vel = site_init.vel + 0.5 * dt * (site.acc - site_init.vel * chi)
                site.angvel = site_init.angvel + 0.5 * dt * (site.angacc - site_init.angvel * chi)
            if abs(chi - chi_old) < EPSILON:
                break
        if i == n_max:
            print "NVT UPDATE DID NOT CONVERGE\n\n"

        self.data['nvt_chi'] = chi
        self.data['nvt_chi_dt'] += 0.5 * dt * chi
        return
        
        
    def update_step(self):
        """ select different integrator """
        integrator = self.config.give(keyword="integrator")
        if integrator == "NVE":
            self.update_step_nve()
            
        return
        
    def print_info(self):
        """
        print dynamics info.
        """
        pot = self.status.give(keyword="potential_energy")
        kin = self.get_kinetic_energy()
        ene = self.get_energy_nve()
        temperature = self.get_temperature()        
        print("%30s %16.10lf\n" % ("POTENTIAL ENERGY", pot))
        print("%30s %16.10lf\n" % ("KINETIC ENERGY", kin))
        print("%30s %16.10lf\n" % ("INVARIANT", ene))
        print("%30s %16.10lf\n" % ("TEMPERATURE", temperature))
        
        return
        
    def print_restart(self, fp):
        """ restart dump """
        print "RESTART DATA DUMP"
        sites = self.status.dump(fp)  
        
        
        return

     
        
    def print_status(self, fp):
        print("    INITIAL STATE\n\n")
        self.print_restart(fp)
        self.print_info()
        return
        
    def kill(self):
        # finalize...
        
        return

    def sim(self):
        """ md worker """
        print("MOLECULAR DYNAMICS JOB\n\n\n")
        fp = open("trans.xyz", "w")
        self.setup()
        self.velocitize()
        self.remove_system_drift()
        self.action()
        self.print_status(fp)
        
        n_max_steps = self.config.give(keyword="max_steps")
        n_print_steps = self.config.give(keyword="print_steps")
        for i in xrange(n_max_steps):
            self.update_step()
            if i % n_print_steps == 0:
                self.print_status(fp)                
        self.kill()        
        fp.close()
        print("MOLECULAR DYNAMICS JOB COMPLETED SUCCESSFULLY\n")

        return

    
    
if __name__ == "__main__":
    config = bConfig()
    dyn = bDynSOC(config)
    dyn.setup()
    dyn.dump()
    
    
    # dyn.sim()
    
    # dyn.setup()
    # dyn.velocitize()
    # dyn.get_kinetic_energy()
    # t = dyn.get_temperature()
    # print dyn.get_energy_nve()
    # print dyn.get_system_com()
    # print dyn.get_system_cov()
    # print dyn.get_system_angmom()
    # print dyn.get_system_inertia_tensor()
    # dyn.remove_system_drift()
    # dyn.update_step_nvt()
    # print dyn.get_energy_nve()
    