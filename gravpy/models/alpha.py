from basemodel import BaseModel
from sie import elliptical, spherical


class Alpha(BaseModel):
    def __init__(self, b, x0, y0, e, te, s, alpha):
        super(Alpha, self).__init__(x0, y0, e, te)
        self.b = b
        self.alpha = alpha
        # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.
        self.s = s if s != 0.0 else 1e-4

    def modelargs(self, alpha=False):
        if alpha:
            return [self.b, self.x0, self.y0, self.e, self.te, self.s, self.alpha]
        else:
            return [self.b, self.x0, self.y0, self.e, self.te, self.s]

    @BaseModel.standard_frame_rotation
    def phiarray(self, x, y, numexpr=True, *args, **kwargs):
        modelargs = self.modelargs()
        modelargs_with_alpha = self.modelargs(alpha=True)

        if self.alpha == 1.0:
            if self.e == 0.0:
                return spherical(x, y, modelargs, numexpr=numexpr)
            else:
                return elliptical(x, y, modelargs, numexpr=numexpr)
        elif self.alpha == -1.0:
            return plummer(x, y, modelargs)
        else:
            try:
                from fortran.alpha import alphaf
                return alphaf.general(x, y, modelargs_with_alpha)
            except:
                raise Exception("Alpha!=(0 | -1) not implemented yet")


def plummer(x, y, modelargs):
    '''Calculation for the alpha=-1 case'''
    b, x0, y0, e, te, s = modelargs

    x2 = x * x
    y2 = y * y
    s2 = s * s
    q = 1.0 - e
    q2 = q * q
    om = 1.0 - q2
    rt = np.sqrt(om)
    front = b ** 3 * q / s
    psi = np.sqrt(q2 * (s2 + x2) + y2)
    psi2 = psi * psi
    psi3 = psi2 * psi

    invDenom = 1 / (psi * (om + x2 + (psi + s) ** 2))
    phix = front * x * (psi + s * q2) * invDenom
    phiy = front * y * (psi + s) * invDenom

    tmp = (psi + s) ** 2 + om * x2

    # phixx
    phixx = -2.0 * front * x2 * ((psi + q2 * s) / tmp) ** 2 / psi2 + front * (
    psi2 * (psi + q2 * s) - q2 ** 2 * x2 * s) / (psi3 * tmp)
    phiyy = -2.0 * front * y2 * ((psi + s) / tmp) ** 2 / psi2 + front * (psi2 * (psi + s) - y2 * s) / (psi3 * tmp)
    phixy = -front * x * y * (2.0 * (psi + s) * (psi + q2 * s) / (psi2 * tmp * tmp) + q2 * s / (psi3 * tmp))
    pot = 0.5 * front * np.log((psi * s) ** 2 + om * x2) - front * np.log(s * (1.0 + np.fabs(q)))

    return np.array((pot, phix, phiy, phixx, phiyy, phixy))


def general(x, y, modelargs):
    b, x0, y0, e, te, s, a = modelargs

    def single_eval(x, y):
        r = np.sqrt(x * x + y * y)
        cost = x / r if r > 0. else 1.0
        sint = y / r if r > 0. else 0.0
        front = (1 - e) / 2.0 * np.power(b, 2.0 - a)

        if abs(a) > 0.01:
            front *= a

        if abs(e) < 1.0e-6:
            phir = front * alpha0phir(r, s, a)
            phir_r = front * np.power(s, a - 2.0) if r == 0.0 else phir / r
            phirr = -phir_r + 2.0 * front * np.power(s * s + r * r, 0.5 * a - 1.0)
            phix = phir_r * x
            phiy = phir_r * y
            phixx = phir_r * sint * sint + phirr * cost * cost
            phiyy = phir_r * cost * cost + phirr * sint * sint
            phixy = (phirr - phir_r) * sint * cost

            if abs(s) <= 2.0 * r:
                warnings.warn("Calling hypergeo and gamma functions")
                pot = sp.hyp2f1(-0.5 * a, -0.5 * a, 1.0 - 0.5 * a, -(s * s) / (r * r))
                pot *= b * b * np.power(r / b, a) / (a * a)
                pot -= b * b * np.power(s / b, a) / a * (np.log(r / s) + 0.5 * (0.577216 - sp.digamma(-0.5 * a)))
                if abs(a) > 0.01:
                    pot *= a

            else:
                pot = 0.0 if r == 0.0 else front * scipyint.quad(alpha0phir, 0.0, r, args=(s, a))

        else:
            p0, p1, p2, p3, p4, p5 = alpha_integral(0.0, 1.0, r, e, sint, cost, s, a)
            pot = front * (p0 / 2.)
            phix = front * (p1 * x)
            phiy = front * (p2 * y)
            phixx = front * (p1 + 2. * x * x * p3)
            phiyy = front * (p2 + 2. * y * y * p4)
            phixy = front * (2. * x * y * p5)

        return pot, phix, phiy, phixx, phiyy, phixy

    return np.transpose([single_eval(x, y) for x, y in itertools.izip(x, y)])


def alpha0phir(r, s, a):
    if r / s < 1.0e-4:
        return r * np.power(s, a - 2.0)
    elif a == 0.0:
        return np.log(1.0 + r * r / (s * s)) / r
    else:
        return 2.0 / (a * r) * (np.power(s * s + r * r, a / 2.0) - np.power(s, a))


def alpha_integral(lower, upper, r, e, sint, cost, s, a):
    def function(u, num):
        q = 1.0 - e
        t0 = 1.0 - (1.0 - q * q) * u
        t1 = 1.0 / np.sqrt(t0)
        t3 = t1 / t0
        t5 = t3 / t0
        t6 = cost * cost + sint * sint / (1.0 - (1.0 - q * q) * u) * r * r
        t7 = t6 * u

        mphiu = t6 * np.power(s, a - 2.0) if u == 0.0 else np.sqrt(t7) * alpha0phir(np.sqrt(t7), s, a) / u
        k = np.power(s * s + t7, a / 2.0 - 1.0)
        kp = k * (a / 2.0 - 1.0) / (s * s + t7)

        # mphiu,k,kp = alphakap(u,e,sint,cost,s,a)

        if num == 0:
            return t1 * mphiu
        elif num == 1:
            return t1 * k
        elif num == 2:
            return t3 * k
        elif num == 3:
            return t1 * kp * u
        elif num == 4:
            return t5 * kp * u
        elif num == 5:
            return t3 * kp * u

    return [scipyint.quad(function, lower, upper, args=(num))[0] for num in range(6)]

