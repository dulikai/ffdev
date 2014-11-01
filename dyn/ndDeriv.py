#! python 

class ndDeriv:
    """ numerical differical """
    def modify():
        self.cluster['site'][isite]['coord'][idim] += delta
        return
    def update():
        
        return
    def nd(self, delta):
        delta = 1.0e-6
        n_site = self.cluster['n_site']
        site = self.cluster['site']
        for i_site in xrange(n_site):
            mass = site[i_site]['mass']
            acc = np.array([0.0, 0.0, 0.0])
            for j_dim in xrange(3):
                self.modify_site_coord(i_site, j_dim, delta)
                v2 = self.potential(sigma, epsilon)
                self.modify_site_coord(i_site, j_dim, -2*delta)
                v1 = self.potential(sigma, epsilon)
                self.modify_site_coord(i_site, j_dim, delta)
                acc[j_dim] = - (v2 - v1) / (2.0 * delta) / mass
            site[i_site]['acc'] = acc
            # print "acc",acc
