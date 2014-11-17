#! /usr/bin/env python

# information control simulations.

class bConfig():
    """
    configuration 
    """
    time_step = 0.5        # fs
    temperature = 298.15    # K
    enable_pbc = "no"       # pbc condition; not implemented
    ld_gamma = 100          # suppose to be ps^-1
    ld_temperature = 300    # K [default: temperature]
    thermostat_tau = 1.0    # unit
    #
    max_steps = 1000
    print_steps = 1
    integrator = "NVE"
    #
    quantum = {
                'n_state': 2,
                'i_state': 2,
                'n_spin': 2,
                'i_spin': 3,
                'qm_method': 'DFT'
              }
              
              
    def __init__(self):
    
        return
		
    def give(self, keyword=""):
        try:
            value = getattr(self, keyword)
        except AttributeError:
            print("@error: no such configuration keyword!!!")
            exit(1)
        return value
        
        
        

if __name__ == "__main__":
    bc = bConfig()
    print bConfig.time_step
    print bc.give(keyword = "enable_pbc")
    
    