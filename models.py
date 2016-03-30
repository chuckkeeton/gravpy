from abc import abstractmethod, ABCMeta
import pysie, pyalpha
from math import sin, cos, pi
import numpy as np
from functools import wraps
from integration import Integrator
from collections import Iterable


class BaseModel(object):
    __metaclass__ = ABCMeta

    def __init__(self, x0, y0, e, te, s):
        """
        Arguments:
            x0 {number} -- x position
            y0 {number} -- y position
            e {number} -- ellipticity
            te {number} -- [description]
            s {number} -- [description]
        """
        self.x0 = x0
        self.y0 = y0
        self.e = e
        self.te = te
        # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.
        self.s = s if s != 0.0 else 1e-4

    @abstractmethod
    def phiarray(self, x, y, *args, **kwargs):
        """Method that returns an array of
        (phi, dphi/dx, dphi/dy, d^2phi/dx^2, d^2phi/dy^2, d^2phi/(dxdy))
        at the given coordinates x,y, where phi is the gravitational potential
        """
        pass


def standard_frame_rotation(phiarray_function):
    """A wrapper that rotates the incoming x,y values into the standard
    frame for lensing calculations and rotates the phiarray values back
    into the frame they were orginally in.
    """

    @wraps(phiarray_function)
    def rotation(self, x, y, *args, **kwargs):
        x0, y0 = self.x0, self.y0
        te = self.te

        c = cos(te * pi / 180)
        c2 = c * c
        s = sin(te * pi / 180)
        s2 = s * s
        sc = s * c

        # transformation from actual coords to natural coords
        # (the frame/axes are rotated so there is no ellipticity angle in the calculation).
        # This makes the expressions in the model modules simpler to calculate.
        xp = -s * (x - x0) + c * (y - y0)
        yp = -c * (x - x0) - s * (y - y0)

        pot, px, py, pxx, pyy, pxy = phiarray_function(self, xp, yp, *args, **kwargs)

        # Inverse transformation back into desired coordinates.
        new_phix = -s * px - c * py
        new_phiy = c * px - s * py
        new_phixx = s2 * pxx + c2 * pyy + 2 * sc * pxy
        new_phiyy = c2 * pxx + s2 * pyy - 2 * sc * pxy
        new_phixy = sc * (pyy - pxx) + (s2 - c2) * pxy

        return np.array((pot, new_phix, new_phiy, new_phixx, new_phiyy, new_phixy))

    return rotation


class SIE(BaseModel):
    def __init__(self, b, x0, y0, e, te, s):
        super(SIE, self).__init__(x0, y0, e, te, s)
        self.b = b

    @standard_frame_rotation
    def phiarray(self, x, y, numexpr=True, *args, **kwargs):
        modelargs = [self.b, self.x0, self.y0, self.e, self.te, self.s]

        if self.e == 0:
            return pysie.spherical(x, y, modelargs, numexpr=numexpr)
        else:
            return pysie.elliptical(x, y, modelargs, numexpr=numexpr)


class Alpha(BaseModel):
    def __init__(self, b, x0, y0, e, te, s, alpha):
        super(Alpha, self).__init__(x0, y0, e, te, s)
        self.b = b
        self.alpha = alpha

    @standard_frame_rotation
    def phiarray(self, x, y, numexpr=True, *args, **kwargs):
        modelargs = [self.b, self.x0, self.y0, self.e, self.te, self.s]

        if self.alpha == 1.0:
            if self.e == 0.0:
                return pysie.spherical(x, y, modelargs, numexpr=numexpr)
            else:
                return pysie.elliptical(x, y, modelargs, numexpr=numexpr)
        elif self.alpha == -1.0:
            return pyalpha.plummer(x, y, modelargs)
        else:
            raise Exception("Alpha!=(0 | -1) not implemented yet")


class NFW(BaseModel):
    """Navarro-Frenk-White profile"""

    def __init__(self, b, x0, y0, e, te, s, ks, rs):
        super(NFW, self).__init__(x0, y0, e, te, s)
        self.b = b
        self.ks = ks
        self.rs = rs
        self.q = 1.0 - self.e
        self.integrator = Integrator(self.q, self.kappa, self.kappa_prime, self.phi_r)

    @standard_frame_rotation
    def phiarray(self, x, y, *args, **kwargs):
        if not isinstance(x, Iterable) and not isinstance(y, Iterable):
            x = [x]
            y = [y]

        potential, phix, phiy, phixx, phiyy, phixy = [[] for i in xrange(6)]
        for i, (local_x, local_y) in enumerate(zip(x, y)):
            print "on {0} out of {1}".format(i, len(x))
            potential.append(0)
            phix.append(self.integrator.phi_x(local_x, local_y))
            phiy.append(self.integrator.phi_y(local_x, local_y))
            phixx.append(self.integrator.phi_xx(local_x, local_y))
            phiyy.append(self.integrator.phi_yy(local_x, local_y))
            phixy.append(self.integrator.phi_xy(local_x, local_y))
        return np.array((potential, phix, phiy, phixx, phiyy, phixy))

    @staticmethod
    def funcF(x):
        dx = x**2 - 1.0
        if abs(x) < 1e-2:
            # series with O(x^6) error
            log2x = np.log(2.0 / x)
            return log2x + x**2 * (0.5 * log2x - 0.25) * x**4 * (0.375 * log2x - 0.21875)
        elif abs(dx) < 1e-2:
            # series with O(dx^6) error
            return 1.0 - (dx / 3.0) + (dx**2 / 5.0) - (dx**3 / 7.0) + (dx**4 / 9.0) - (dx**5 / 11.0)
        elif x > 1.0:
            tmp = np.sqrt(x**2 - 1.0)
            return np.arctan(tmp) / tmp
        else:
            tmp = np.sqrt(1.0 - x**2)
            return np.arctanh(tmp) / tmp

    @staticmethod
    def funcF_prime(x):
        return (1.0 - x**2 * NFW.funcF(x)) / (x * (x**2 - 1.0))

    def kappa(self, r):
        x = r / self.rs
        return 2.0 * self.ks * ((1.0 - NFW.funcF(x)) / (x ** 2.0 - 1.0))

    def kappa_prime(self, r):
        x = r / self.rs
        numerator = (2.0 * self.rs * self.ks * ((self.rs**2 - r**2) * NFW.funcF_prime(x)
                                                + 2.0 * self.rs * r * NFW.funcF(x)
                                                - 2.0 * self.rs * r))
        denomenator = (self.rs**2 - r**2)**2
        return numerator / denomenator

    def phi(self, r):
        x = r / self.rs
        return 2.0 * self.ks * self.rs ** 2 * (np.log(x / 2.0) ** 2 - np.arctanh(np.sqrt(1.0 - x ** 2)) ** 2)

    def phi_r(self, x):
        """Spherical deflection"""
        return 4.0 * self.ks * self.rs * ((np.log(x / 2.0) + self.funcF(x)) / x)
