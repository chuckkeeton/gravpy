import abc
from math import sin,cos,pi
import numpy as np
from functools import wraps

from .sie.python import siepy
from .alpha.python import alphapy
# f2py modules imported below, clunky pathing due to fortran modules, full function path is fmodel.fmodel.routine
from .sie.fortran.sief import sief
from .alpha.fortran.alphaf import alphaf 
from .nfw.fortran.nfwf import nfwf

class baseModel:
    __metaclass__ = abc.ABCMeta
    @abc.abstractmethod
    def modelargs(self):
        '''Define and return a list of the model parameters'''
        
    @abc.abstractmethod
    def phiarray(self,x,y):
        '''Method that returns an array of [phi, dphi/dx, dphi/dy, d^2phi/dx^2, d^2phi/dy^2, d^2phi/(dxdy)] at the given coordinates x,y, where phi is the gravitational potential'''

def standard_frame_rotation(phiarray_function):
    '''A wrapper that rotates the incoming x,y values into the standard frame for lensing calculations and rotates the phiarray values back into the frame they were orginally in.'''
    @wraps(phiarray_function)
    def rotation(self,x,y,*args,**kwargs):
        x0,y0 = self.x0,self.y0
        te = self.te
        
        c = cos(te*pi/180);c2=c*c
        s = sin(te*pi/180);s2=s*s;sc=s*c
        
        #transformation from actual coords to natural coords (the frame/axes are rotated so there is no ellipticity angle in the calculation). This makes the expressions in the model modules simpler to calculate.
        xp = -s*(x-x0) + c*(y-y0)
        yp = -c*(x-x0) - s*(y-y0)

        pot,px,py,pxx,pyy,pxy = phiarray_function(self,xp,yp,*args,**kwargs) #get the phiarray values
        
        # Inverse transformation back into desired coordinates.
        new_phix = -s*px-c*py
        new_phiy = c*px-s*py 
        new_phixx= s2*pxx+c2*pyy+2*sc*pxy
        new_phiyy= c2*pxx+s2*pyy-2*sc*pxy
        new_phixy= sc*(pyy-pxx)+(s2-c2)*pxy
        
        return np.array((pot,new_phix,new_phiy,new_phixx,new_phiyy,new_phixy))
    
    return rotation
            

class SIE(baseModel):
    
    def __init__(self,b,x0,y0,e,te,s):
        for key,value in locals().items():
            setattr(self, key, value)
        self.s = s if s!=0.0 else 1e-4 # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.
    def modelargs(self):
        return [self.b,self.x0,self.y0,self.e,self.te,self.s]
    
    @standard_frame_rotation
    def phiarray(self,x,y,numexpr=True):
        modelargs = self.modelargs()
        
        if self.e==0:
            return siepy.spherical(x,y,modelargs,numexpr=numexpr)
        else:
            return siepy.elliptical(x,y,modelargs,numexpr=numexpr)


class alpha(baseModel):
    def __init__(self,b,x0,y0,e,te,s,alpha):
        for key,value in locals().items():
            setattr(self, key, value)
        self.s = s if s!=0.0 else 1e-4 # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.

    def modelargs(self,alpha=False):
        if alpha:
            return [self.b,self.x0,self.y0,self.e,self.te,self.s,self.alpha]
        else:
            return [self.b,self.x0,self.y0,self.e,self.te,self.s]
        
    @standard_frame_rotation
    def phiarray(self,x,y,numexpr=True):
        modelargs = self.modelargs()
        modelargs_with_alpha = self.modelargs(alpha=True)
        
        if self.alpha==1.0:
            if self.e == 0.0:
                return siepy.spherical(x,y,modelargs,numexpr=numexpr)
            else:
                return siepy.elliptical(x,y,modelargs,numexpr=numexpr)
        elif self.alpha== -1.0:
            return alphapy.plummer(x,y,modelargs)
        else:
            return alphaf.general(x,y,modelargs_with_alpha)

class nfw(baseModel):
    def __init__(self,b,x0,y0,e,te,rs):
        for key,value in locals().items():
            setattr(self, key, value)
        if self.rs == 0.0:
            raise ValueError("NFW scale radius cannot be zero")

    def modelargs(self):
        return [self.b,self.x0,self.y0,self.e,self.te,self.rs]

    @standard_frame_rotation
    def phiarray(self,x,y,numexpr=True):
        modelargs = self.modelargs()
        
        return nfwf.nfw(x,y,modelargs)

        

    
