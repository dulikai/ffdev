
struct nvt_data {
	double chi;
	double chi_dt;
};
    
    def update_step_nvt(self):
        """
        update md step in nvt...
        the Nose-Hoover thermostat is used.
        ref: 
        Canonical dynamics: Equilibrium phase-space distributions;
        William G. Hoover, Phys. Rev. A 31, 1695 (1985)
        """
        dt = self.get_md_config(keyword = "time_step")
        target = self.get_md_config(keyword = "temperature")
        tau = self.get_md_config(keyword = "thermostat_tau")
        nvt = self.get_md_config(keword = "nvt_cfg")
        
        t0 = self.get_temperature()
        n_site = self.get_md_system(keyword = "n_site")
        cluster = self.get_md_system(keyword = "cluster")
        
        for site in cluster:
            site.vel += 0.5 * dt * (site.acc - site.vel * nvt.chi)
            site.pos += site.vel * dt
        
        # reset nvt_chi and nvt_chi_dt
        nvt.chi += 0.5 * dt * (t0 / target - 1.0) / tau /tau
        nvt.chi_dt += 0.5 * dt * nvt.chi
        ffmodel.compute_forces(cluster)
        
        chi_init = nvt.chi
        vel_init = []
        for site in cluster:
            vel_init.append(site.vel)
        
        max_iter = self.get_md_config(keyword = "max_iter")
        for iter in xrange(max_iter):
            chi_old = nvt.chi
            t = self.get_temperature()
            ratio = t / targe
            # update chi
            nvt.chi = chi_init + 0.5 * dt * (ratio - 1.0) / tau / tau
            
            for site in cluster:
                site.vel = vel_init + 0.5 * dt * (site.acc - vel_init * nvt.chi)
            
            if abs(nvt.chi - chi_old) < EPSILON):
                break
                
            if iter == max_iter:
                print "WARNING: NVT UPDATE DID NOT CONVERGE\n\n"
                
	nvt->chi_dt += 0.5 * dt * nvt->chi;
    return
    
    def get_kinetic_energy(self):
        """ summ ke """
        ke = 0.0
        cluster = self.get_md_system(keyword = "cluster")
        
        for site in cluster:
            ke += site.mass * np.dot(site.vel, site.vel)
            
        return 0.5 * ke
    

    def get_temperature(self):
        """ T system """
        ke = self.get_kinetic_energy()
        t = 2.0 * ke / BOLTZMANN / n_freedom
        return
        

    def update_step_npt(self):
        """
        update md step in nvt...
        the Hoover thermostat is used.
        ref: 
        Hoover NPT dynamics for systems varying in shape and size
        Simone Melchionna, Giovanni Ciccotti, Brad Lee Holian
        Mol. Phys. 78, 533 (1993)
        """
        npt = get_md_config(keyword = "npt_cfg")
        
        dt = self.get_md_config(keyword = "time_step")
        t_target = self.get_md_config(keyword = "temperature")
        t_tau = self.get_md_config(keyword = "thermostat_tau")
        p_target = self.get_md_config(keyword = "pressure")
        p_tau = self.get_md_config(keyword = "barostat_tau")

        t_tau2 = t_tau * t_tau
        p_tau2 = p_tau * p_tau
        kbt = BOLTZMANN * t_target
        
	double t_tau2 = t_tau * t_tau;
	double p_tau2 = p_tau * p_tau;
	double kbt = BOLTZMANN * t_target;

	double t0 = get_temperature(md);
	double p0 = get_pressure(md);
	double v0 = get_volume(md);

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * dt * (body->force.x / body->mass -
					body->vel.x * (data->chi + data->eta));
		body->vel.y += 0.5 * dt * (body->force.y / body->mass -
					body->vel.y * (data->chi + data->eta));
		body->vel.z += 0.5 * dt * (body->force.z / body->mass -
					body->vel.z * (data->chi + data->eta));

		body->angmom.x += 0.5 * dt * (body->torque.x -
						body->angmom.x * data->chi);
		body->angmom.y += 0.5 * dt * (body->torque.y -
						body->angmom.y * data->chi);
		body->angmom.z += 0.5 * dt * (body->torque.z -
						body->angmom.z * data->chi);

		rotate_body(body, dt);
	}

        return
        
 struct npt_data {
	double chi;
	double chi_dt;
	double eta;
};

static void update_step_npt(struct md *md)
{
	struct npt_data *data = (struct npt_data *)md->data;

	double dt = cfg_get_double(md->cfg, "time_step");
	double t_tau = cfg_get_double(md->cfg, "thermostat_tau");
	double t_target = cfg_get_double(md->cfg, "temperature");
	double p_tau = cfg_get_double(md->cfg, "barostat_tau");
	double p_target = cfg_get_double(md->cfg, "pressure");

	double t_tau2 = t_tau * t_tau;
	double p_tau2 = p_tau * p_tau;
	double kbt = BOLTZMANN * t_target;

	double t0 = get_temperature(md);
	double p0 = get_pressure(md);
	double v0 = get_volume(md);

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * dt * (body->force.x / body->mass -
					body->vel.x * (data->chi + data->eta));
		body->vel.y += 0.5 * dt * (body->force.y / body->mass -
					body->vel.y * (data->chi + data->eta));
		body->vel.z += 0.5 * dt * (body->force.z / body->mass -
					body->vel.z * (data->chi + data->eta));

		body->angmom.x += 0.5 * dt * (body->torque.x -
						body->angmom.x * data->chi);
		body->angmom.y += 0.5 * dt * (body->torque.y -
						body->angmom.y * data->chi);
		body->angmom.z += 0.5 * dt * (body->torque.z -
						body->angmom.z * data->chi);

		rotate_body(body, dt);
	}

	data->chi += 0.5 * dt * (t0 / t_target - 1.0) / t_tau2;
	data->chi_dt += 0.5 * dt * data->chi;
	data->eta += 0.5 * dt * v0 * (p0 - p_target) / md->n_bodies / kbt / p_tau2;

	vec_t com = get_system_com(md);
	vec_t pos_init[md->n_bodies];

	for (int i = 0; i < md->n_bodies; i++)
		pos_init[i] = md->bodies[i].pos;

	for (int iter = 1; iter <= MAX_ITER; iter++) {
		bool done = true;

		for (int i = 0; i < md->n_bodies; i++) {
			struct body *body = md->bodies + i;
			vec_t pos = wrap(md, &body->pos);

			vec_t v = {
				data->eta * (pos.x - com.x),
				data->eta * (pos.y - com.y),
				data->eta * (pos.z - com.z)
			};

			vec_t new_pos = {
				pos_init[i].x + dt * (body->vel.x + v.x),
				pos_init[i].y + dt * (body->vel.y + v.y),
				pos_init[i].z + dt * (body->vel.z + v.z)
			};

			done = done && vec_dist(&body->pos, &new_pos) < EPSILON;
			body->pos = new_pos;
		}

		if (done)
			break;

		if (iter == MAX_ITER)
			printf("WARNING: NPT UPDATE DID NOT CONVERGE\n\n");
	}

	vec_scale(&md->box, exp(dt * data->eta));
	check_fail(efp_set_periodic_box(md->efp, md->box.x, md->box.y, md->box.z));

	compute_forces(md);

	double chi_init = data->chi, eta_init = data->eta;
	vec_t angmom_init[md->n_bodies], vel_init[md->n_bodies];

	for (int i = 0; i < md->n_bodies; i++) {
		angmom_init[i] = md->bodies[i].angmom;
		vel_init[i] = md->bodies[i].vel;
	}

	for (int iter = 1; iter <= MAX_ITER; iter++) {
		double chi_prev = data->chi;
		double eta_prev = data->eta;
		double t_cur = get_temperature(md);
		double p_cur = get_pressure(md);
		double v_cur = get_volume(md);

		data->chi = chi_init + 0.5 * dt * (t_cur / t_target - 1.0) / t_tau2;
		data->eta = eta_init + 0.5 * dt * v_cur * (p_cur - p_target) /
							md->n_bodies / kbt / p_tau2;

		for (int i = 0; i < md->n_bodies; i++) {
			struct body *body = md->bodies + i;

			body->vel.x = vel_init[i].x + 0.5 * dt *
				(body->force.x / body->mass -
					vel_init[i].x * (data->chi + data->eta));
			body->vel.y = vel_init[i].y + 0.5 * dt *
				(body->force.y / body->mass -
					vel_init[i].y * (data->chi + data->eta));
			body->vel.z = vel_init[i].z + 0.5 * dt *
				(body->force.z / body->mass -
					vel_init[i].z * (data->chi + data->eta));

			body->angmom.x = angmom_init[i].x + 0.5 * dt *
				(body->torque.x - angmom_init[i].x * data->chi);
			body->angmom.y = angmom_init[i].y + 0.5 * dt *
				(body->torque.y - angmom_init[i].y * data->chi);
			body->angmom.z = angmom_init[i].z + 0.5 * dt *
				(body->torque.z - angmom_init[i].z * data->chi);
		}

		if (fabs(data->chi - chi_prev) < EPSILON &&
		    fabs(data->eta - eta_prev) < EPSILON)
			break;

		if (iter == MAX_ITER)
			printf("WARNING: NPT UPDATE DID NOT CONVERGE\n\n");
	}

	data->chi_dt += 0.5 * dt * data->chi;
}
            