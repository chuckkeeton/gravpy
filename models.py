from abc import abstractmethod, ABCMeta
import pysie, pyalpha
from math import sin, cos, pi
import numpy as np
from functools import wraps


class BaseModel(object):
    __metaclass__ = ABCMeta

    def __init__(self, x0, y0, e, te, s):
        """[summary]
        
        [description]
        
        Arguments:
            x0 {float} -- x position
            y0 {float} -- y position
            e {float} -- ellipticity
            te {float} -- [description]
            s {float} -- [description]
        """
        self.x0 = x0
        self.y0 = y0
        self.e = e
        self.te = te
        # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.
        self.s = s if s != 0.0 else 1e-4

    @abstractmethod
    def phiarray(self, x, y):
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
    def phiarray(self, x, y, numexpr=True):
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
    def phiarray(self, x, y, numexpr=True):
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
