"""Real-kernel Hartley transforms, with the unnormalized book convention."""
import numpy as np


def hartley(x, axis=-1):
    """Complex-linear H, valid for real or complex data of any length."""
    x = np.asarray(x)
    F = np.fft.fft(x, axis=axis)
    reverse = (-np.arange(x.shape[axis])) % x.shape[axis]
    H = ((1+1j)*F + (1-1j)*np.take(F, reverse, axis=axis))/2
    return H.real if np.isrealobj(x) else H


def hartley_twiddle(x, shift=.5):
    x = np.asarray(x)
    n = x.size
    angle = 2*np.pi*shift*np.arange(n)/n
    return np.cos(angle)*x + np.sin(angle)*x[(-np.arange(n)) % n]


def fht(x):
    """Radix-two fast Hartley transform; no complex arithmetic for real x."""
    x = np.asarray(x)
    n = x.size
    if x.ndim != 1 or n == 0 or n & (n-1):
        raise ValueError('Use a nonempty vector of power-of-two length.')
    if n == 1:
        return x.copy()
    even = fht(x[::2])
    odd = hartley_twiddle(fht(x[1::2]))
    return np.concatenate((even+odd, even-odd))


def hartley_convolution(a, b):
    """Circular convolution of equal-length vectors using Hartley transforms."""
    a, b = np.asarray(a), np.asarray(b)
    if a.ndim != 1 or a.shape != b.shape or a.size == 0:
        raise ValueError('The inputs must be nonempty vectors of equal length.')
    c, d = hartley(a), hartley(b)
    reverse = (-np.arange(a.size)) % a.size
    spectrum = (c*(d+d[reverse]) + c[reverse]*(d-d[reverse]))/2
    return hartley(spectrum)/a.size


def fractional_hartley(x, power):
    """A fixed-branch unitary power of H/sqrt(N), not of unnormalized H."""
    x = np.asarray(x)
    transformed = hartley(x)/np.sqrt(x.size)
    return (x+transformed)/2 + np.exp(1j*np.pi*power)*(x-transformed)/2


def hartley2d(x):
    x = np.asarray(x)
    if x.ndim != 2:
        raise ValueError('Use a two-dimensional array.')
    return hartley(hartley(x, axis=0), axis=1)


def finite_fht(x, root, iota):
    """Radix-two H over an odd finite field; root has order len(x), iota^2=-1."""
    field, n = type(x), x.size
    if n == 0 or n & (n-1):
        raise ValueError('The length must be a power of two.')
    if root**n != 1 or (n > 1 and root**(n//2) == 1) or iota**2 != -field(1):
        raise ValueError('The root orders or square root of -1 are incorrect.')
    if n == 1:
        return x.copy()
    even, odd = finite_fht(x[::2], root**2, iota), finite_fht(x[1::2], root**2, iota)
    a = root**np.arange(n//2)
    cosine, sine = (a+a**-1)/field(2), (a-a**-1)/(field(2)*iota)
    mixed = cosine*odd + sine*odd[(-np.arange(n//2)) % (n//2)]
    return field(np.concatenate((even+mixed, even-mixed)))


def finite_hartley_example():
    import galois
    prime = galois.GF(5, compile='python-calculate')
    polynomial = galois.Poly([1, 0, 0, 0, 3], field=prime)
    field = galois.GF(5**4, irreducible_poly=polynomial, compile='python-calculate')
    x = field([0, 1, 1, 1, 2, 2, 1, 4, 0, 1, 2, 0, 1, 0, 2, 0])
    root, iota = field(5), field(2)  # X modulo X^4+3, and a square root of -1
    y = finite_fht(x, root, iota)
    assert np.array_equal(finite_fht(y, root, iota), x)
    return x, y
