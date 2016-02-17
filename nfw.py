import numpy as np

def density(ps, r, rs):
	return ps / ((r / rs) * (1.0 + (r / rs)))**2.0

def kappa(r, Ks, x):
	return  2.0 * Ks * ((1.0 - F(x)) / (x**2.0 - 1.0))

def F(x):
	if x == 1.0:
		return 1
	else if x > 1.0:
		return (1.0 / np.sqrt(x**2 - 1.0)) * np.arctan(np.sqrt(x**2 - 1.0))
	else:
		return (1.0 / np.sqrt(1.0 - x**2)) * np.arctanh(1.0 - x**2)

def phi(Ks, rs, x):
	return 2.0 * Ks * rs**2 * (np.log(x / 2.0)**2 - np.arctanh(np.sqrt(1.0 - x**2))**2)

def phi_r():
	return 4.0 * Ks * rs * ((np.log(x / 2.0) + F(x)) / x)