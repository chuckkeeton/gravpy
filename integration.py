import scipy.integrate as integrate

def f7(x, y):
    q = 0.5
    j0 = jn(0)
    result = j0(x, y, q)
    return ((q * x * result[0]), result[1])

def f8(x, y):
    q = 0.5
    j1 = jn(1)
    result = j1(x, y, q)
    return ((q * x * result[0]), result[1])

def kappa(e):
    alpha = 1.0
    b = 2.0
    s = 0.01
    return (0.5 * b**(2 - alpha)) / ((s**2 + e**2)**(1 - (alpha / 2.0)))

def xi_squared(u, x, y, q):
    return u * (x**2 + (y**2 / (1.0 - ((1.0  - q**2) * u))))

def jn(n):
    def local_jn(u, x, y, q):
        return (kappa(xi_squared(u, x, y, q)) / ((1.0 - (1.0 - q**2) * u)**(n + 0.5)))
    # return lambda x, y, q: integrate.quad(lambda u: local_jn(u, x, y, q), 0, 1)  # ~550 micro-seconds
    return lambda x, y, q: integrate.quad(local_jn, 0, 1, args=(x, y, q))  # ~510-520 micro-seconds
    #return lambda x, y, q: integrate.quadrature(local_jn, 0, 1, args=(x, y, q))  # ~3.2 - 3.4 ms
    #return lambda x, y, q: integrate.romberg(local_jn, 0, 1, args=(x, y, q))  # ~9.35 - 9.38 ms

print f7(1, 1)
print f8(1, 1)