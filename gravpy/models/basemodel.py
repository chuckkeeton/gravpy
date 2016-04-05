from abc import ABCMeta, abstractmethod
from functools import wraps
import numpy as np


class BaseModel(object):
    __metaclass__ = ABCMeta

    def __init__(self, x0, y0, e, te):
        """
        Arguments:
            x0 {number} -- x position
            y0 {number} -- y position
            e {number} -- ellipticity
            te {number} -- theta_ellipticity
        """
        self.x0 = x0
        self.y0 = y0
        self.e = e
        self.te = te

    @abstractmethod
    def phiarray(self, x, y, *args, **kwargs):
        """Method that returns an array of
        (phi, dphi/dx, dphi/dy, d^2phi/dx^2, d^2phi/dy^2, d^2phi/(dxdy))
        at the given coordinates x,y, where phi is the gravitational potential
        """
        pass

    @staticmethod
    def standard_frame_rotation(phiarray_function):
        """A wrapper that rotates the incoming x,y values into the standard
        frame for lensing calculations and rotates the phiarray values back
        into the frame they were orginally in.
        """

        @wraps(phiarray_function)
        def rotation(self, x, y, *args, **kwargs):
            x0, y0 = self.x0, self.y0
            te = self.te

            c = np.cos(te * np.pi / 180)
            c2 = c * c
            s = np.sin(te * np.pi / 180)
            s2 = s * s
            sc = s * c

            # transformation from actual coords to natural coords
            # (the frame/axes are rotated so there is no ellipticity angle in the calculation).
            # This makes the expressions in the model modules simpler to calculate.
            xp = -s * (x - x0) + c * (y - y0)
            yp = -c * (x - x0) - s * (y - y0)

            pot, px, py, pxx, pyy, pxy = phiarray_function(self, xp, yp, *args, **kwargs)  # get the phiarray values

            # Inverse transformation back into desired coordinates.
            new_phix = -s * px - c * py
            new_phiy = c * px - s * py
            new_phixx = s2 * pxx + c2 * pyy + 2 * sc * pxy
            new_phiyy = c2 * pxx + s2 * pyy - 2 * sc * pxy
            new_phixy = sc * (pyy - pxx) + (s2 - c2) * pxy

            return np.array((pot, new_phix, new_phiy, new_phixx, new_phiyy, new_phixy))

        return rotation
