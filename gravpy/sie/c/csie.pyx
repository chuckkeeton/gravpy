# cython: profile=True
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt,log,atan,atanh


def phiarray(xi,yi,modelargs,vec=False):
    if vec:
        return phiarray_vec(xi,yi,modelargs)
    else:
        return phiarray_for(xi,yi,modelargs)


cpdef phiarray_for(np.ndarray[np.float64_t, ndim=1] xi, 
                   np.ndarray[np.float64_t, ndim=1] yi, 
                   np.ndarray[np.float64_t, ndim=1] modelargs):

    #xi.shape ~ (n), modelargs.shape ~ (parameter_array_length)
    #output.shape ~ (6,n)
    cdef Py_ssize_t xsize = xi.shape[0]

    cdef np.ndarray[np.float64_t, ndim=2] output = np.empty((xsize,6),dtype=np.float64)
        
    which_function(xi,yi,modelargs,output,xsize)
    
    return np.transpose(output)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef phiarray_vec(np.ndarray[np.float64_t, ndim=2] xi, 
                   np.ndarray[np.float64_t, ndim=2] yi, 
                   np.ndarray[np.float64_t, ndim=3] modelargs):

    #xi.shape ~ (n,num_of_mass_components), modelargs.shape ~ (parameter_array_length,1,num_of_mass_components)
    #output.shape ~ (6,n,num_of_mass_components)
    cdef Py_ssize_t xsize = xi.shape[0]
    cdef Py_ssize_t num_masscomp = modelargs.shape[2]
    cdef Py_ssize_t param_len = modelargs.shape[0]
        
        
    cdef np.ndarray[np.float64_t, ndim=2] new_modelargs = np.squeeze(modelargs,axis=1).T.copy() #shape~ (num_masscomp,param_length)
    cdef np.ndarray[np.float64_t, ndim=2] new_x = xi.T.copy() #(num_of_mass_components,n)
    cdef np.ndarray[np.float64_t, ndim=2] new_y = yi.T.copy() #(num_of_mass_components,n)
    cdef np.ndarray[np.float64_t, ndim=3] output = np.empty((num_masscomp,xsize,6),dtype=np.float64)
        
    cdef unsigned int i
    
    for i in range(num_masscomp):
        which_function(new_x[i],new_y[i],new_modelargs[i],output[i],xsize)
            
            
    return np.transpose(output,[2,1,0])

    

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void which_function(double[:] x, double[:] y, double[:] model_line, double[:,:] output, unsigned int xsize):
    cdef double e = model_line[3]
    
    cdef unsigned int j
    if e == 0.0:
        for j in range(xsize):
            spherical(x[j],y[j],model_line,output[j])
    else:
        for j in range(xsize):
            elliptical(x[j],y[j],model_line,output[j])

    #(x,y,model_line)

    
@cython.cdivision(True)    
cdef void elliptical(double x, double y, double[:] modelargs, double[:] output):
    cdef double b  = modelargs[0]
    cdef double x0 = modelargs[1]
    cdef double y0 = modelargs[2]
    cdef double e  = modelargs[3]
    cdef double te = modelargs[4]
    cdef double s  = modelargs[5]
    if s == 0.0:
        s = 1e-4

    cdef double x2  = x**2
    cdef double y2  = y**2
    cdef double s2  = s**2
    cdef double q   = 1.0-e
    cdef double q2  = q**2
    cdef double om  = 1.0-q2
    cdef double rt  = sqrt(om)
    cdef double psi = sqrt(q2*(s2+x2)+y2)
    cdef double psis= psi + s

    cdef double phix = b*q/rt *atan(rt*x/psis)
    cdef double phiy = b*q/rt *atanh(rt*y/(psi+s*q2))

    cdef double invDenom = 1/(psi*(om*x2+psis**2))
    cdef double phixx = b*q*(psi*psis-q2*x2)*invDenom
    cdef double phiyy = b*q*(x2+s*psis)*invDenom
    cdef double phixy = -b*q*x*y*invDenom

    cdef double pot = b*q*s*(-0.5*log(psis**2+om*x2) + log(s*(1.0+q)) ) + x*phix+y*phiy
    
    output[0] = pot
    output[1] = phix
    output[2] = phiy
    output[3] = phixx
    output[4] = phiyy
    output[5] = phixy



@cython.cdivision(True)    
cdef void spherical(double x, double y, double[:] modelargs, double[:] output):
    cdef double b  = modelargs[0]
    cdef double x0 = modelargs[1]
    cdef double y0 = modelargs[2]
    cdef double e  = modelargs[3]
    cdef double te = modelargs[4]
    cdef double s  = modelargs[5]
    if s == 0.0:
        s = 1e-4
    
    cdef double r = sqrt(x**2+y**2)
    cdef double rad = sqrt(r**2+s**2)
    cdef double sprad = s + rad
    cdef double invDenom = 1/(rad*sprad**2)
    
    cdef double pot = b * (rad-s*(1+log(sprad/(2*s))))
    cdef double phix = b * x / sprad
    cdef double phiy = b * y / sprad
    cdef double phixx = b * (s*sprad + y**2) * invDenom
    cdef double phiyy = b * (s*sprad + x**2) * invDenom
    cdef double phixy = -b*x*y *invDenom
    
    output[0] = pot
    output[1] = phix
    output[2] = phiy
    output[3] = phixx
    output[4] = phiyy
    output[5] = phixy


