#! /usr/bin/env python

#
# config site
struct md {
	size_t n_bodies;
	struct body *bodies;
	size_t n_freedom;
	vec_t box;
	double potential_energy;
	double (*get_invariant)(const struct md *);
	void (*update_step)(struct md *);
	struct state *state;
	void *data;
};


class DynInfo:
    """
    dynamics simulation data
    config, site, pbc, state
    """
    def __init__(self):
        self.site = bSite()
        self.n_sites = 0
        self.box = []
        self.config = bConfig()        
        return
        
        
        
        
        
if __name__ == "__main__":
    info = DynInfo()
    