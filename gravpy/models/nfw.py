import numpy as np
from basemodel import BaseModel
from collections import Iterable
from integration import Integrator


class NFW(BaseModel):
    """Navarro-Frenk-White profile"""

    def __init__(self, b, x0, y0, e, te, s, ks, rs):
        super(NFW, self).__init__(x0, y0, e, te)
        self.b = b
        self.ks = ks
        self.rs = rs
        self.q = 1.0 - self.e
        self.integrator = Integrator(self.q, self.kappa, self.kappa_prime, self.phi_r)
        # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.
        self.s = s if s != 0.0 else 1e-4

    @BaseModel.standard_frame_rotation
    def phiarray(self, x, y, *args, **kwargs):
        if not isinstance(x, Iterable) and not isinstance(y, Iterable):
            x = [x]
            y = [y]

        potential, phix, phiy, phixx, phiyy, phixy = [[] for i in xrange(6)]
        print "calculating phiarray for {0} (x, y) pairs".format(len(x))
        for i, (local_x, local_y) in enumerate(zip(x, y)):
            potential.append(0)
            phix.append(self.integrator.phi_x(local_x, local_y))
            phiy.append(self.integrator.phi_y(local_x, local_y))
            phixx.append(self.integrator.phi_xx(local_x, local_y))
            phiyy.append(self.integrator.phi_yy(local_x, local_y))
            phixy.append(self.integrator.phi_xy(local_x, local_y))
        return np.array((potential, phix, phiy, phixx, phiyy, phixy))

    @staticmethod
    def funcF(x):
        dx = x ** 2 - 1.0
        if abs(x) < 1e-2:
            # series with O(x^6) error
            log2x = np.log(2.0 / x)
            return log2x + x ** 2 * (0.5 * log2x - 0.25) * x ** 4 * (0.375 * log2x - 0.21875)
        elif abs(dx) < 1e-2:
            # series with O(dx^6) error
            return 1.0 - (dx / 3.0) + (dx ** 2 / 5.0) - (dx ** 3 / 7.0) + (dx ** 4 / 9.0) - (dx ** 5 / 11.0)
        elif x > 1.0:
            tmp = np.sqrt(x ** 2 - 1.0)
            return np.arctan(tmp) / tmp
        else:
            tmp = np.sqrt(1.0 - x ** 2)
            return np.arctanh(tmp) / tmp

    @staticmethod
    def funcF_prime(x):
        return (1.0 - x ** 2 * NFW.funcF(x)) / (x * (x ** 2 - 1.0))

    def kappa(self, r):
        x = r / self.rs
        return 2.0 * self.ks * ((1.0 - NFW.funcF(x)) / (x ** 2.0 - 1.0))

    def kappa_prime(self, r):
        x = r / self.rs
        numerator = (2.0 * self.rs * self.ks * ((self.rs ** 2 - r ** 2) * NFW.funcF_prime(x)
                                                + 2.0 * self.rs * r * NFW.funcF(x)
                                                - 2.0 * self.rs * r))
        denomenator = (self.rs ** 2 - r ** 2) ** 2
        return numerator / denomenator

    def phi(self, r):
        x = r / self.rs
        return 2.0 * self.ks * self.rs ** 2 * (np.log(x / 2.0) ** 2 - np.arctanh(np.sqrt(1.0 - x ** 2)) ** 2)

    def phi_r(self, x):
        """Spherical deflection"""
        return 4.0 * self.ks * self.rs * ((np.log(x / 2.0) + self.funcF(x)) / x)
