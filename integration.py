from scipy.integrate import quad
from numpy import sqrt


def phi(x, y, q, phi_r):
    """Lensing potential"""
    local_i = i(phi_r, q)
    result, err = local_i(x, y)
    return (q / 2.0) * result


def phi_x(x, y, q, kappa):
    """deflection in x at position (x, y)
    
    Arguments:
        x {number} -- x coordinate
        y {number} -- y coordinate
        q {number} -- q = 1 - elipticity
        kappa {function} -- takes in 1 argument: xi
                            (eliptical coordinate)
    
    Returns:
        number -- the deflection in x
    """
    j0 = jn(0, q, kappa)
    result, err = j0(x, y)
    return q * x * result


def phi_y(x, y, q, kappa):
    j1 = jn(1, q, kappa)
    result, err = j1(x, y)
    return q * x * result


def phi_xx(x, y, q, kappa, kappa_prime):
    k0 = kn(0, q, kappa_prime)
    j0 = jn(0, q, kappa)
    result_k, err_k = k0(x, y)
    result_j, err_j = j0(x, y)
    return 2.0 * q * x ** 2.0 * result_k + q * result_j


def phi_yy(x, y, q, kappa, kappa_prime):
    k2 = kn(2, q, kappa_prime)
    j1 = jn(1, q, kappa)
    result_k, err_k = k2(x, y)
    result_j, err_j = j1(x, y)
    return 2.0 * q * y ** 2.0 * result_k + q * result_j


def phi_xy(x, y, q, kappa_prime):
    k1 = kn(1, q, kappa_prime)
    result, err = k1(x, y)
    return 2.0 * q * x * y * result


def sple_kappa(w, alpha, b, s):
    # w = xi**2,
    # so xi**2 is not necessary
    return (0.5 * b ** (2 - alpha)) / ((s ** 2 + w) ** (1 - (alpha / 2.0)))


def sple_kappa_prime(w, alpha, b, s):
    # the derivative of kappa over w (xi**2)
    return 0.25 * (alpha - 2.0) * b ** (2.0 - alpha) * (s ** 2.0 + w) ** ((alpha / 2.0) - 2.0)


def xi(x, y, q):
    return sqrt(xi_squared(x, y, q))


def xi_squared(x, y, q):
    return x**2 + (y**2 / q**2)


def xi_u(u, x, y, q):
    return sqrt(xi_u_squared(u, x, y, q))


def xi_u_squared(u, x, y, q):
    return u * (x ** 2 + (y ** 2 / (1.0 - ((1.0 - q ** 2) * u))))


def i(phi_r, q):
    def integrand(u, x, y):
        return (xi_u(u, x, y, q) / u) * phi_r(xi(u, x, y, q)) / (1.0 - (1.0 - q**2) * u)**0.5
    return lambda x, y: quad(integrand, 0, 1, args=(x, y))


def jn(n, q, kappa):
    def integrand(u, x, y):
        return kappa(xi_u_squared(u, x, y, q)) / ((1.0 - (1.0 - q**2) * u)**(n + 0.5))
    return lambda x, y: quad(integrand, 0, 1, args=(x, y))


def kn(n, q, kappa_prime):
    def integrand(u, x, y):
        return u * kappa_prime(xi_u_squared(u, x, y, q)) / ((1.0 - (1.0 - q**2) * u)**(n + 0.5))
    return lambda x, y: quad(integrand, 0, 1, args=(x, y))

    ###########
    # analytic_answers = sie.elliptical(1, 1, [2.0, None, None, 0.5, None, 0.01])

    # print phi_x(1, 1)
    # print phi_y(1, 1)

    # print phi_xx(1, 1)
    # print analytic_answers[3]

    # print phi_yy(1, 1)
    # print analytic_answers[4]

    # print phi_xy(1, 1)
    # print analytic_answers[5]
    ###########
    #
