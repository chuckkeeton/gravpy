import numpy as np
import numexpr as ne
import core
#modelargs: (major) radius, x-center position, y-center position, ellipticity, ellipticity angle, core radius

def phiarray(xi,yi,modelargs,vec=False,numexpr=True):
    np.place(modelargs[-1],modelargs[-1]==0,1e-4) # replaces core radius (s)==0 -> 1e-4, fixes /0 situations in potential calculation. 
    b,x0,y0,e,te,s  = modelargs[:6]
    
    if vec:
        return core.cond_break(xi,yi,modelargs,[e==0,e!=0],[spherical,elliptical])
    else:
        
        if e==0:
            return spherical(xi,yi,modelargs,numexpr=numexpr)
        else:
            return elliptical(xi,yi,modelargs,numexpr=numexpr)
    
        
def elliptical(x,y,modelargs,numexpr=True):
    b,x0,y0,e,te,s  = modelargs[:6]
    
    x2  = x*x
    y2  = y*y
    s2  = s*s
    q   = 1.0-e
    q2  = q*q
    om  = 1.0-q2
    rt  = np.sqrt(om)
    psi = np.sqrt(q2*(s2+x2)+y2) if not numexpr else ne.evaluate("sqrt(q2*(s2+x2)+y2)")
    psis= psi + s

    if numexpr:
        phix = ne.evaluate("b*q/rt *arctan(rt*x/psis)")
        phiy = ne.evaluate("b*q/rt *arctanh(rt*y/(psi+s*q2))")

        invDenom = ne.evaluate("1/(psi*(om*x2+psis*psis))")
        phixx = ne.evaluate("b*q*(psi*psis-q2*x2)*invDenom")
        phiyy = ne.evaluate("b*q*(x2+s*psis)*invDenom")
        phixy = ne.evaluate("-b*q*x*y*invDenom")

        pot = ne.evaluate("b*q*s*(-0.5*log(psis*psis+om*x2) + log(s*(1.0+q)) ) + x*phix+y*phiy")

    else:
        phix = b*q/rt *np.arctan(rt*x/psis)
        phiy = b*q/rt *np.arctanh(rt*y/(psi+s*q2))
        
        invDenom = 1/(psi*(om*x2+psis*psis))
        phixx = b*q*(psi*psis-q2*x2)*invDenom
        phiyy = b*q*(x2+s*psis)*invDenom
        phixy = -b*q*x*y*invDenom

        pot = b*q*s*(-0.5*np.log(psis*psis+om*x2) + np.log(s*(1.0+q)) ) + x*phix+y*phiy
    
    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

def spherical(x,y,modelargs,numexpr=True):
    b,x0,y0,e,te,s  = modelargs[:6]
    
    if numexpr:
        rad = ne.evaluate("sqrt(x*x+y*y+s*s)")
        sprad = ne.evaluate("s + rad")
        invDenom = ne.evaluate("1/(rad*sprad*sprad)")
    
        pot = ne.evaluate("b * (rad-s*(1+log(sprad/(2*s))))")
        phix = ne.evaluate("b * x / sprad")
        phiy = ne.evaluate("b * y / sprad")
        phixx = ne.evaluate("b * (s*sprad + y*y) * invDenom")
        phiyy = ne.evaluate("b * (s*sprad + x*x) * invDenom")
        phixy = ne.evaluate("-b*x*y *invDenom")
    
    else:
        rad = np.sqrt(x*x+y*y+s*s)
        sprad = s + rad
        invDenom = 1/(rad*sprad*sprad)
    
        pot = b * (rad-s*(1+np.log(sprad/(2*s))))
        phix = b * x / sprad
        phiy = b * y / sprad
        phixx = b * (s*sprad + y*y) * invDenom
        phiyy = b * (s*sprad + x*x) * invDenom
        phixy = -b*x*y *invDenom

    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

