import numpy as np
import numexpr as ne

# modelargs: (major) radius, x-center position, y-center position, ellipticity, ellipticity angle, core radius


if ne.use_vml:
    ne.set_vml_accuracy_mode('fast')


def elliptical(x, y, modelargs, numexpr=True):
    b, x0, y0, e, te, s = modelargs[:6]

    x2 = x * x
    y2 = y * y
    s2 = s * s
    q = 1.0 - e
    q2 = q * q
    om = 1.0 - q2
    rt = np.sqrt(om)
    psi = np.sqrt(q2 * (s2 + x2) + y2) if not numexpr else ne.evaluate("sqrt(q2*(s2+x2)+y2)")
    psis = psi + s

    if numexpr:
        phix = ne.evaluate("b*q/rt *arctan(rt*x/psis)")
        phiy = ne.evaluate("b*q/rt *arctanh(rt*y/(psi+s*q2))")

        invDenom = ne.evaluate("b*q/(psi*(om*x2+psis*psis))")
        phixx = ne.evaluate("(psi*psis-q2*x2)*invDenom")
        phiyy = ne.evaluate("(x2+s*psis)*invDenom")
        phixy = ne.evaluate("-x*y*invDenom")

        pot = ne.evaluate("b*q*s*(-0.5*log(psis*psis+om*x2) + log(s*(1.0+q)) ) + x*phix+y*phiy")

    else:
        phix = b * q / rt * np.arctan(rt * x / psis)
        phiy = b * q / rt * np.arctanh(rt * y / (psi + s * q2))

        invDenom = b * q / (psi * (om * x2 + psis * psis))
        phixx = (psi * psis - q2 * x2) * invDenom
        phiyy = (x2 + s * psis) * invDenom
        phixy = -x * y * invDenom

        pot = b * q * s * (-0.5 * np.log(psis * psis + om * x2) + np.log(s * (1.0 + q))) + x * phix + y * phiy

    return np.array((pot, phix, phiy, phixx, phiyy, phixy))


def spherical(x, y, modelargs, numexpr=True):
    b, x0, y0, e, te, s = modelargs[:6]

    if numexpr:
        rad = ne.evaluate("sqrt(x*x+y*y+s*s)")
        sprad = ne.evaluate("s + rad")
        invDenom = ne.evaluate("b/(rad*sprad*sprad)")

        pot = ne.evaluate("b * (rad-s*(1+log(sprad/(2*s))))")
        phix = ne.evaluate("b * x / sprad")
        phiy = ne.evaluate("b * y / sprad")
        phixx = ne.evaluate("(s*sprad + y*y) * invDenom")
        phiyy = ne.evaluate("(s*sprad + x*x) * invDenom")
        phixy = ne.evaluate("-x*y *invDenom")

    else:
        rad = np.sqrt(x * x + y * y + s * s)
        sprad = s + rad
        invDenom = b / (rad * sprad * sprad)

        pot = b * (rad - s * (1 + np.log(sprad / (2 * s))))
        phix = b * x / sprad
        phiy = b * y / sprad
        phixx = (s * sprad + y * y) * invDenom
        phiyy = (s * sprad + x * x) * invDenom
        phixy = -x * y * invDenom

    return np.array((pot, phix, phiy, phixx, phiyy, phixy))
