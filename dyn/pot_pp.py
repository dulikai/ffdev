#! python

# mono- point charge potential 

class pot_pp():
    def __init__(self):
        return
    def energy(self, r, Q1, Q2):
        """
        V_elst = Q_1 Q_2 / r
        """
        s = Q1 * Q2 / r
        s6 = s**6; s12 = s6 * s6
        pot = 4.0 * epsilon * (s12 - s6)
        return pot
        
    def gradient(self):
        """
        V' = 4 \epsilon [-12/r x^12 - (-6/r) x^6]
        """
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        grad = 4.0 * epsilon * ((-12.0/r) * s12 - (-6/r) * s6)
        return grad
    
    def eandg(self, ri, rj, rk, sigma, epsilon):
        """
        energy & gradient
        """
        r = np.linalg.norm(ri-rj)
        s = sigma / r
        s6 = s**6; s12 = s6 * s6
        ene = 4.0 * epsilon * (s12 - s6)
        grad = 4.0 * epsilon * ((-12.0/r) * s12 - (-6/r) * s6)
        return ene, grad
        
        