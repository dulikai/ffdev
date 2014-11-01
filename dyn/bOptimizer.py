#! /usr/bin/env python

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html


from scipy.optimize import minimize

# optimizer

class bOptimizer():
    """
    multivariable minimizer
    """
    
    def __init__(self):
        
        return
        
    def cons_modify(self):
        """
        modify the force by constraints
        """
        bCons.give()
        
        return
    def wash(self):
        """
        optimizer
        """
        options = {'gtol': 1e-6, 'disp': True}
        res = minimize(pes, vars, grad, options)
        

        return
        