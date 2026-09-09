"""Centered splines and cardinal kernels in the book's Fourier convention."""
import numpy as np
from scipy.interpolate import BSpline


def centered_bspline(degree, x):
    x = np.asarray(x)
    if degree < 0 or int(degree) != degree:
        raise ValueError('Spline degree must be a nonnegative integer.')
    if degree == 0:
        return ((x >= -.5) & (x < .5)).astype(float)
    knots = np.arange(degree+2)-(degree+1)/2
    values = BSpline.basis_element(knots, extrapolate=False)(x)
    return np.nan_to_num(values, nan=0)


def spline_symbol(degree, omega):
    radius = int(np.ceil((degree+1)/2))
    k = np.arange(-radius, radius+1)
    return np.cos(np.asarray(omega)[..., None]*k) @ centered_bspline(degree, k)


def cardinal_spline_spectrum(degree, omega):
    omega = np.asarray(omega)
    return np.sinc(omega/(2*np.pi))**(degree+1)/spline_symbol(degree, omega)


def cardinal_cubic(x):
    x = np.asarray(x)
    k = np.floor(x)[..., None].astype(int)+np.arange(-1, 3)
    coefficients = np.sqrt(3)*(np.sqrt(3)-2)**abs(k)
    return np.sum(coefficients*centered_bspline(3, x[..., None]-k), axis=-1)


def mirror_initial(u, q):
    """Exact causal initial value for whole-sample mirror extension."""
    u = np.asarray(u)
    n = u.size
    if n == 1:
        return u[0]/(1-q)
    j = np.arange(1, n-1)
    numerator = u[0] + q**(n-1)*u[-1] + np.sum((q**j+q**(2*n-2-j))*u[j])
    return numerator/(1-q**(2*n-2))


def cubic_prefilter_parallel(u):
    """Cubic B-spline coefficients, using two parallel mirror-boundary recursions."""
    u = np.asarray(u)
    if u.ndim != 1 or u.size == 0:
        raise ValueError('Use a nonempty vector.')
    q = np.sqrt(3)-2
    forward = np.empty(u.size, dtype=np.result_type(u, float))
    backward = np.empty_like(forward)
    forward[0], backward[-1] = mirror_initial(u, q), mirror_initial(u[::-1], q)
    for k in range(1, u.size):
        forward[k] = u[k] + q*forward[k-1]
    for k in range(u.size-2, -1, -1):
        backward[k] = u[k] + q*backward[k+1]
    return np.sqrt(3)*(forward+backward-u)


def cubic_prefilter_cascade(u):
    """The same mirror-boundary prefilter as a causal/anticausal cascade."""
    u = np.asarray(u)
    if u.ndim != 1 or u.size == 0:
        raise ValueError('Use a nonempty vector.')
    if u.size == 1:
        return u.astype(np.result_type(u, float))
    q = np.sqrt(3)-2
    d = np.empty(u.size, dtype=np.result_type(u, float))
    d[0] = mirror_initial(u, q)
    for k in range(1, u.size):
        d[k] = u[k] + q*d[k-1]
    c = np.empty_like(d)
    c[-1] = -6*q*(d[-1]+q*d[-2])/(1-q*q)
    for k in range(u.size-2, -1, -1):
        c[k] = -6*q*d[k] + q*c[k+1]
    return c
