import numpy as np

def plummer(x,y,modelargs):
    '''Calculation for the alpha=-1 case'''
    b,x0,y0,e,te,s  = modelargs
    
    x2 = x*x
    y2 = y*y
    s2 = s*s
    q  = 1.0-e
    q2 = q*q
    om = 1.0-q2
    rt = np.sqrt(om)
    front = b**3*q/s
    psi = np.sqrt(q2*(s2+x2)+y2)
    psi2 = psi*psi
    psi3 = psi2*psi
    
    invDenom = 1/(psi*(om+x2+(psi+s)**2))
    phix = front*x*(psi+s*q2)*invDenom
    phiy = front*y*(psi+s)   *invDenom
    
    tmp = (psi+s)**2+om*x2

    # phixx
    phixx = -2.0*front*x2*((psi+q2*s)/tmp)**2/psi2 + front*(psi2 *(psi+q2*s)-q2**2*x2*s)/(psi3*tmp)
    
    phiyy = -2.0*front*y2*((psi+   s)/tmp)**2/psi2 + front*(psi2 *(psi   +s)-      y2*s)/(psi3*tmp)

    phixy = -front*x*y*(2.0*(psi+s)*(psi+q2*s)/(psi2*tmp*tmp) + q2*s/(psi3*tmp))

    pot = 0.5*front*np.log((psi*s)**2+om*x2) - front*np.log(s*(1.0+np.fabs(q)))
    
    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

    
    

