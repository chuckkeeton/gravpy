import numpy as np

#modelargs: (major) radius, x-center position, y-center position, ellipticity, ellipticity angle, core radius

def phiarray(xi,yi,modelargs):
    np.place(modelargs[-1],modelargs[-1]==0,0.0001) # replaces core radius (s)==0 -> 0.0001, fixes /0 situations in pot calculation. 
    b,x0,y0,e,te,s  = modelargs

    # following is some filtering and ordering to seperate elliptical cases from spherical cases and then recombine resulting phiarrays back into the same order they came in argument-wise
    #TODO: make wrapper for generalizing this block of code, expect this kind of condition break to be useful in additional model routines...
    n = np.max(xi.shape)
    empty_shape = np.zeros((6,n,0),dtype='float64')

    where_e0     = np.flatnonzero(e==0) 
    where_e_not0 = np.flatnonzero(e!=0) 
    
    sphericalargs = np.take(modelargs,where_e0,axis=2)
    spheremodels = spherical(xi[:,where_e0],yi[:,where_e0],sphericalargs) if sphericalargs.size!= 0 else empty_shape
    ellipticalargs = np.take(modelargs,where_e_not0,axis=2)
    ellipticalmodels = elliptical(xi[:,where_e_not0],yi[:,where_e_not0],ellipticalargs) if ellipticalargs.size!=0 else empty_shape
         
    allmodels = np.concatenate((spheremodels,ellipticalmodels),axis=2)
    sorted_indices = np.argsort(np.hstack((where_e0,where_e_not0)))
    sorted_models = np.take(allmodels,sorted_indices,axis=2)

    return sorted_models
    
        
def elliptical(x,y,modelargs):
    b,x0,y0,e,te,s  = modelargs
    
    x2  = x**2
    y2  = y**2
    s2  = s**2
    q   = 1.0-e
    q2  = q**2
    om  = 1.0-q2
    rt  = np.sqrt(om)
    psi = np.sqrt(q2*(s2+x2)+y2)
    psis= psi + s

    phix = b*q/rt *np.arctan(rt*x/psis)
    phiy = b*q/rt *np.arctanh(rt*y/(psi+s*q2))

    invDenom = 1/(psi*(om*x2+psis**2))
    phixx = b*q*(psi*psis-q2*x2)*invDenom
    phiyy = b*q*(x2+s*psis)*invDenom
    phixy = -b*q*x*y*invDenom

    pot = b*q*s*(-0.5*np.log(psis**2+om*x2) + np.log(s*(1.0+q)) ) + x*phix+y*phiy


    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

def spherical(x,y,modelargs):
    b,x0,y0,e,te,s  = modelargs

    r = np.sqrt(x**2+y**2)
    rad = np.sqrt(r**2+s**2)
    sprad = s + rad
    invDenom = 1/(rad*sprad**2)
    

    pot = b * (rad-s*(1+np.log(sprad/(2*s))))
    phix = b * x / sprad
    phiy = b * y / sprad
    phixx = b * (s*sprad + y**2) * invDenom
    phiyy = b * (s*sprad + x**2) * invDenom
    phixy = -b*x*y *invDenom

    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

