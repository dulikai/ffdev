#! /usr/bin/env python

import os
import numpy as np

import copy

class QuantumInt:
    """ the quantum amplitue propagator """
    def __init__(self):
        h_old = np.array([[1001.0, 0.0], [0.0, 1001.9]])
        h_new = np.array([[1000.2, 0.0], [0.0, 1002.8]])
        h_ref = np.array([[1000., 0.], [0.0, 1000.0]])
        rho = np.array([[1., 0.], [0.0, 0.0]])
        nac_time_old = np.array([[0., 0.8], [1.2, 0.0]])
        nac_time_new = np.array([[0., 0.7], [1.0, 0.0]])
        
        self.vars = {
            "n_elec": 100,
            "n_state": 2,
            "dtime": 0.01,
            "h_old": h_old,
            "h_new": h_new,
            "h_ref": h_ref,
            "nac_type": "time",
            "nac_time_old": nac_time_old,
            "nac_time_new": nac_time_new,
            "rho": rho
            }
        
        return
        
    # f' = f(t)
    # f(t) = 0.5 * t
    # f(t+\delta t) = f(t) + f'(t+0.5*\delta t) \delta t
    def test(self):
        n = 1000
        dt = 0.001
        t = 0.0
        fold = 0.0
        fgrad = 0.5 * (dt/2)# * (dt/2)
        fnew = fold + fgrad * dt
        for i in xrange(n):
            t = t + i*dt
            fold = fnew
            fgrad = 0.5 * (t+dt/2)#*(t+dt/2)
            fnew = fold + fgrad * dt
            print t, fnew
        return

        
        
    @staticmethod    
    def func(x, y):
        f = -1 * x * y * y
        return f
        
    def rkf4(self):
        # input a,b,x0,y0,n:0 5 0 2 20
        a = 0.
        b = 5.
        x0 = 0.
        y0 = 2.
        n = 20
        print "x0\ty0\tk1\tk2\tk3\tk4\n"
        h = (b - a) / n
        for i in xrange(n):
            k1 = self.func(x0, y0)
            k2 = self.func(x0+h/2, y0+k1*h/2)
            k3 = self.func(x0+h/2, y0+k2*h/2)
            k4 = self.func(x0+h, y0+h*k3)
            # print "%12.6f%12.6f" % (x0, y0),
            # print "%12.6f%12.6f" % (k1, k2),
            # print "%12.6f%12.6f\n" % (k3, k4)
            y0 += h*(k1+2*k2+2*k3+k4)/6
            x0 += h
        print "EULER 4: %12.6f%12.6f" % (x0, y0)    


    def rkf1(self):
        a = 0.
        b = 5.
        x0 = 0.
        y0 = 2.
        n = 20
        print "x0\ty0\tk1\tk2\tk3\tk4\n"
        h = (b - a) / n
        for i in xrange(n):
            k1 = self.func(x0, y0)
            # print "%12.6f%12.6f" % (x0, y0),
            # print "%12.6f" % (k1)
            y0 += h*(k1)
            x0 += h
        print "EULER 1: %12.6f%12.6f" % (x0, y0)    
        return
    
    def rkf2(self):
        """
        yi+1=yi+h*(K1+K2)/2
        K1=f(xi,yi)  K2=f(xi+h,yi+h*K1)
        """
        a = 0.
        b = 5.
        x0 = 0.
        y0 = 2.
        n = 20
        print "x0\ty0\tk1\tk2\tk3\tk4\n"
        h = (b - a) / n
        for i in xrange(n):
            k1 = self.func(x0, y0)
            k2 = self.func(x0+h, y0+k1*h)
            # print "%12.6f%12.6f" % (x0, y0),
            # print "%12.6f%12.6f" % (k1, k2)
            y0 += h*(k1+k2)/2
            x0 += h
        print "EULER 2: %12.6f%12.6f" % (x0, y0)    
        return
        
    def rkf3(self):
        """
        yi+1=yi+h*(K1+K2)/2
        K1=f(xi,yi)  K2=f(xi+h,yi+h*K1)
        """
        a = 0.
        b = 5.
        x0 = 0.
        y0 = 2.
        n = 20
        print "x0\ty0\tk1\tk2\tk3\tk4\n"
        h = (b - a) / n
        for i in xrange(n):
            k1 = self.func(x0, y0)
            k2 = self.func(x0+h/2., y0+k1*h/2.)
            k3 = self.func(x0+h, y0-k1*h+k2*h*2.)
            # print "%12.6f%12.6f" % (x0, y0),
            # print "%12.6f%12.6f" % (k1, k2)
            y0 += h*(k1+4*k2+k3)/6.0
            x0 += h
        print "EULER 3: %12.6f%12.6f" % (x0, y0)    
        return
         
         
         
         
    #
    # i dc_j/dt = \sum_i {c_i [H_{ji}-\epsilon_0 \delta_{ji}-i v d_{ji}]}
    #
    def build_nac_eff(self, ratio):
        """ build nac """
        nac_type = self.vars['nac_type']
        n_elec = self.vars['n_elec']
        # time derivative
        if nac_type == "time":
            nac_time_old = self.vars['nac_time_old']
            nac_time_new = self.vars['nac_time_new']
            nac_time_eff = nac_time_old + (nac_time_new-nac_time_old) * ratio
        elif nac_type == "cart":
            nac_cart_old = self.vars['nac_cart_old']
            nac_cart_new = self.vars['nac_cart_new']
            vel_old = self.vars['vel_old']
            vel_new = self.vars['vel_new']
            vel_eff = vel_old + (vel_new - vel_old) * ratio
            nac_cart_eff = nac_cart_old + (nac_cart_new - nac_cart_old) * ratio
            nac_time_eff = np.dot(vel_eff, nac_cart_eff)
        # print nac_time_eff    
        return nac_time_eff
            
    def build_hamilton_eff(self, ratio): 
        """ buidl effective hamilton based on linear interpolation """
        n_elec = self.vars['n_elec']
        H_old = self.vars['h_old']
        H_new = self.vars['h_new']
        H_ref = self.vars['h_ref']
        # calc. electronic hamilton, 
        H_eff = H_old + (H_new-H_old) * ratio
        # subtract reference for numerical convenience.
        H_eff = H_eff - H_ref
        # calc Nonadiabatic coupling with respect to time
        nac_time_eff = self.build_nac_eff(ratio)
        # 
        H_eff = H_eff - nac_time_eff
        # multiply -j
        H_eff = -1.0j * H_eff
        # print H_eff
        return H_eff
        

    def unitary_propagator(self):
        """
        Unitary propagator
        suppose rho = (c*) (c)
        rho(t+dt) = exp(-iAdt) rho(t) exp(iAdt)
        """        
        rho = self.vars['rho']
        n_elec = self.vars['n_elec']
        n_state = self.vars['n_state']
        dtime = self.vars['dtime']
        delta = dtime / float(n_elec)
        for i_elec in xrange(n_elec):
            ratio = float(i_elec) / float(n_elec)
            h_eff = self.build_hamilton_eff(ratio)
            w, v = np.linalg.eigh(h_eff)
            vh = np.transpose(v).conj()
            # exp ( -i A dt )
            eval = np.exp(-1.0j * w * delta)
            emat = np.diag(eval) 
            t = np.dot(v, emat)
            e_left = np.dot(t, vh) 
            # exp ( i A dt )
            eval = np.exp(1.0j * w * delta)
            emat = np.diag(eval)
            t = np.dot(v, emat)
            e_right = np.dot(t, vh)            
            # propagate rho
            t = np.dot(e_left, rho)
            rho = np.dot(t, e_right)
            #
            print rho[0][0].real, rho[1][1].real
        self.vars['rho'] = rho
        print "rho: ", rho
        return 
            
         
        # Amat = np.array([[3, 2+1j], [2-1j, 1]])
        # w, v = np.linalg.eigh(Amat)
        # vh = np.transpose(v).conj()
         
         
         
        
    def get_drho(self, elec_time, rho):
        """
        rho' = -i [rho, H]
        note H*
        """
        dtime = self.vars['dtime']
        n_state = self.vars['n_state']
        ratio = elec_time / float(dtime)
        h_eff = self.build_hamilton_eff(ratio)
        drho = -1.0  * (np.dot(h_eff, rho) - np.dot(rho, h_eff))
        
        # drho = np.zeros((2,2))
        # for k in xrange(n_state):
            # for l in xrange(n_state):
                # a = 0.0
                # for j in xrange(n_state):
                    # a = a + (h_eff[j][l]*rho[k][j] - h_eff[k][j]*rho[j][l])
                # drho[k][l] = - a * 1.0j    
        # print drho
        return drho #np.transpose(drho).conj()
        

        
    def rkf_propagator(self):
        """
        R.K. propagator
        suppose rho = (c*) (c)
        rho(t+dt) = exp(-iAdt) rho(t) exp(iAdt)
        """        
        rho = self.vars['rho']
        
        n_elec = self.vars['n_elec']
        n_state = self.vars['n_state']
        dtime = self.vars['dtime']
        h = dtime / float(n_elec)
        x0 = 0.0
        y0 = copy.deepcopy(rho)
        for i_elec in xrange(n_elec):
            k1 = self.get_drho(x0, y0)
            k2 = self.get_drho(x0+h/2, y0+k1*h/2)
            k3 = self.get_drho(x0+h/2, y0+k2*h/2)
            k4 = self.get_drho(x0+h, y0+h*k3)
            # print "%12.6f%12.6f" % (x0, y0),
            # print "%12.6f%12.6f" % (k1, k2),
            # print "%12.6f%12.6f\n" % (k3, k4)
            # y0 += h*(k1)
            y0 = y0 + h*(k1+2*k2+2*k3+k4)/6
            x0 = x0 + h
            print y0[0][0].real, y0[1][1].real
        

            #
        # print rho[0][0].real, rho[1][1].real
        self.vars['rho'] = y0
        print "rho: ", y0
        return 
     
     

    def get_dc(self, elec_time, c):
        """
        c' = -i A c
        """
        dtime = self.vars['dtime']
        ratio = elec_time / float(dtime)
        h_eff = self.build_hamilton_eff(ratio)
        dc = 1.0 * np.dot(h_eff, c)
        return dc
        
    def rkf_amp(self):
        """
        R.K. propagator
        suppose rho = (c*) (c)
        rho(t+dt) = exp(-iAdt) rho(t) exp(iAdt)
        """        
        n_elec = self.vars['n_elec']
        # n_elec = 1000
        n_state = self.vars['n_state']
        dtime = self.vars['dtime']
        h = dtime / float(n_elec)
        x0 = 0.0
        y0 = np.array([1.0, 0.0])
        for i_elec in xrange(n_elec):
            k1 = self.get_dc(x0, y0)
            k2 = self.get_dc(x0+h/2, y0+k1*h/2)
            k3 = self.get_dc(x0+h/2, y0+k2*h/2)
            k4 = self.get_dc(x0+h, y0+h*k3)
            # print k1
            # print "%12.6f%12.6f" % (x0, y0),
            # print "%12.6f%12.6f" % (k1, k2),
            # print "%12.6f%12.6f\n" % (k3, k4)
            y0 = y0 + h*(k1+2*k2+2*k3+k4)/6
            x0 = x0 + h
            print y0[0].real, y0[1].real
        rho = np.outer(y0, np.transpose(y0).conj())
            #
        print rho[0][0].real, rho[1][1].real
        return 
                        
            
            
# http://zh.wikipedia.org/wiki/%E9%BE%99%E6%A0%BC%EF%BC%8D%E5%BA%93%E5%A1%94%E6%B3%95        
        
# main program
if __name__ == "__main__":
    q = QuantumInt()
    # q.test()
    # q.rkf1()
    # q.rkf2()
    # q.rkf3()
    # q.rkf4()
    q.unitary_propagator()
    # q.rkf_amp()
    # q.rkf_propagator()
    
    