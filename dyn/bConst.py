#! /usr/bin/env python

# http://en.wikipedia.org/wiki/Atomic_units
# http://dynamomd.org/index.php/FAQ



class bConst():
    INFINITY = 1E308
    EPSILON = 1.0e-8
    PI = 3.14159265358979323846
    k_B = 1.0 
    # http://en.wikipedia.org/wiki/Atomic_units#Fundamental_atomic_units
    # TEMPERATURE GIVEN IN ATOMIC UNIT
    au2k = 3.1577464e+5         
    

        
        
if __name__ == "__main__":
    print bConst.au2k