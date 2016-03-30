from functools import wraps
from scipy.integrate import quad
from numpy import sqrt


def memoize(func):
    """Simple decorator to memoize
    
    If running on Python 3.2+, should use
    functools.lru_cache() instead of this.
    """
    results = {}
    @wraps(func)
    def wrapper(*args, **kwargs):
        key = (func.func_name,) + tuple(args) + tuple(kwargs.iteritems())
        if key in results:
            # print "cache hit {0}".format(key)
            return results[key]
        result = func(*args, **kwargs)
        results[key] = result
        return result
    return wrapper


class Integrator(object):
    """Integrator object to solve integrals
    
    This object can be used to solve the integrals for
    phi, phi_x, phi_y, phi_xx, phi_yy, and phi_xy.
    It will automatically cache results from I, J and K integrals
    in order to minimize the number of computations, at the cost
    of space utilized.
    """

    def __init__(self, q, kappa, kappa_prime, phi_r):
        self.q = q
        self.kappa = kappa
        self.kappa_prime = kappa_prime
        self.phi_r = phi_r
    
    def phi(self, x, y):
        """Lensing potential"""
        local_i = self.i(phi_r)
        result, err = local_i(self, x, y)
        return (self.q / 2.0) * result

    def phi_x(self, x, y):
        result, err = self.j0(x, y)
        return self.q * x * result

    def phi_y(self, x, y):
        result, err = self.j1(x, y)
        return self.q * x * result

    def phi_xx(self, x, y):
        result_k, err_k = self.k0(x, y)
        result_j, err_j = self.j0(x, y)
        return 2.0 * self.q * x ** 2.0 * result_k + self.q * result_j

    def phi_yy(self, x, y):
        result_k, err_k = self.k2(x, y)
        result_j, err_j = self.j1(x, y)
        return 2.0 * self.q * y ** 2.0 * result_k + self.q * result_j

    def phi_xy(self, x, y):
        result, err = self.k1(x, y)
        return 2.0 * self.q * x * y * result

    def xi(self, x, y):
        return sqrt(xi_squared(x, y, self.q))

    def xi_squared(self, x, y):
        return x**2 + (y**2 / self.q**2)

    def xi_u(self, u, x, y):
        return sqrt(xi_u_squared(u, x, y, self.q))

    def xi_u_squared(self, u, x, y):
        return u * (x ** 2 + (y ** 2 / (1.0 - ((1.0 - self.q ** 2) * u))))

    def i(self, phi_r):
        def integrand(u, x, y):
            return (self.xi_u(u, x, y) / u) * self.phi_r(self.xi(u, x, y)) / (1.0 - (1.0 - self.q**2) * u)**0.5
        return lambda self, x, y: quad(integrand, 0, 1, args=(x, y))

    @memoize
    def j0(self, x, y):
        local_j0 = self.jn(0)
        return local_j0(x, y)

    @memoize
    def j1(self, x, y):
        local_j1 = self.jn(1)
        return local_j1(x, y)

    @memoize
    def k0(self, x, y):
        local_k0 = self.kn(0)
        return local_k0(x, y)

    @memoize
    def k1(self, x, y):
        local_k1 = self.kn(1)
        return local_k1(x, y)

    @memoize
    def k2(self, x, y):
        local_k2 = self.kn(2)
        return local_k2(x, y)

    def jn(self, n):
        def integrand(u, x, y):
            return self.kappa(self.xi_u_squared(u, x, y)) / ((1.0 - (1.0 - self.q**2) * u)**(n + 0.5))
        return lambda x, y: quad(integrand, 0, 1, args=(x, y))

    def kn(self, n):
        def integrand(u, x, y):
            return u * self.kappa_prime(self.xi_u_squared(u, x, y)) / ((1.0 - (1.0 - self.q**2) * u)**(n + 0.5))
        return lambda x, y: quad(integrand, 0, 1, args=(x, y))

# if __name__ == '__main__':
#     def sple_kappa(w, alpha, b, s):
#         # w = xi**2,
#         # so xi**2 is not necessary
#         return (0.5 * b ** (2 - alpha)) / ((s ** 2 + w) ** (1 - (alpha / 2.0)))


#     def sple_kappa_prime(w, alpha, b, s):
#         # the derivative of kappa over w (xi**2)
#         return 0.25 * (alpha - 2.0) * b ** (2.0 - alpha) * (s ** 2.0 + w) ** ((alpha / 2.0) - 2.0)
    
#     ###########
#     # analytic_answers = sie.elliptical(1, 1, [2.0, None, None, 0.5, None, 0.01])

#     # print phi_x(1, 1)
#     # print phi_y(1, 1)

#     # print phi_xx(1, 1)
#     # print analytic_answers[3]

#     # print phi_yy(1, 1)
#     # print analytic_answers[4]

#     # print phi_xy(1, 1)
#     # print analytic_answers[5]
#     ###########
#     #
