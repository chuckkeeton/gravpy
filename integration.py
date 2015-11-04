import scipy.integrate as integrate
import sie

def phi_x(x, y):
    q = 0.5
    j0 = jn(0)
    result = j0(x, y, q)
    return ((q * x * result[0]), result[1])

def phi_y(x, y):
    q = 0.5
    j1 = jn(1)
    result = j1(x, y, q)
    return ((q * x * result[0]), result[1])

def phi_xx(x, y):
    q = 0.5
    k0 = kn(0)
    j0 = jn(0)
    result_k0 = k0(x, y, q)
    result_j0 = j0(x, y, q)
    return (2.0 * q * x**2.0 * result_k0[0] + q * result_j0[0])

def phi_yy(x, y):
    q = 0.5
    k2 = kn(2)
    j1 = jn(1)
    result_k2 = k2(x, y, q)
    result_j1 = j1(x, y, q)
    return (2.0 * q * y**2.0 * result_k2[0] + q * result_j1[0])

def phi_xy(x, y):
    q = 0.5
    k1 = kn(1)
    result = k1(x, y, q)
    return (2.0 * q * x * y * result[0])

def kappa(w):
    # w = xi**2,
    # so xi**2 is not necessary
    alpha = 1.0
    b = 2.0
    s = 0.01
    return (0.5 * b**(2 - alpha)) / ((s**2 + w)**(1 - (alpha / 2.0)))

def kappa_prime(w):
    # the derivative of kappa over w (xi**2)
    alpha = 1.0
    b = 2.0
    s = 0.01
    return 0.25 * (alpha - 2.0) * b**(2.0 - alpha) * (s**2.0 + w)**((alpha / 2.0) - 2.0)

def xi_squared(u, x, y, q):
    return u * (x**2 + (y**2 / (1.0 - ((1.0  - q**2) * u))))

def jn(n):
    def integrand(u, x, y, q):
        return (kappa(xi_squared(u, x, y, q)) / ((1.0 - (1.0 - q**2) * u)**(n + 0.5)))
    # return lambda x, y, q: integrate.quad(lambda u: integrand(u, x, y, q), 0, 1)  # ~550 micro-seconds
    return lambda x, y, q: integrate.quad(integrand, 0, 1, args=(x, y, q))  # ~510-520 micro-seconds
    #return lambda x, y, q: integrate.quadrature(integrand, 0, 1, args=(x, y, q))  # ~3.2 - 3.4 ms
    #return lambda x, y, q: integrate.romberg(integrand, 0, 1, args=(x, y, q))  # ~9.35 - 9.38 ms

def kn(n):
    def integrand(u, x, y, q):
        return (u * kappa_prime(xi_squared(u, x, y, q)) / ((1.0 - (1.0 - q**2) * u)**(n + 0.5)))
    return lambda x, y, q: integrate.quad(integrand, 0, 1, args=(x, y, q))

analytic_answers = sie.elliptical(1, 1, [2.0, None, None, 0.5, None, 0.01])

print phi_x(1, 1)
print phi_y(1, 1)
print phi_xx(1, 1)
print analytic_answers[3]
print phi_yy(1, 1)
print analytic_answers[4]
print phi_xy(1, 1)
print analytic_answers[5]