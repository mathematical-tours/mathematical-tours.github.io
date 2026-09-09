"""Chirp evaluation, Toeplitz products, and Fourier quadrature."""
import numpy as np


def czt(x, z, count=None):
    """Return sum_j x[j] z**(-j*k), for k=0,...,count-1."""
    x = np.asarray(x, dtype=complex)
    n, m = x.size, x.size if count is None else int(count)
    if x.ndim != 1 or n == 0 or m < 1 or z == 0:
        raise ValueError('Use a nonempty vector, nonzero z, and positive output length.')
    size = 1 << (n+m-2).bit_length()
    logarithm = np.log(complex(z))
    j, k = np.arange(n), np.arange(m)
    differences = np.arange(1-n, m)
    kernel = np.zeros(size, dtype=complex)
    kernel[differences % size] = np.exp(logarithm*differences**2/2)
    h = x*np.exp(-logarithm*j**2/2)
    convolution = np.fft.ifft(np.fft.fft(h, size)*np.fft.fft(kernel))
    return np.exp(-logarithm*k**2/2)*convolution[:m]


def scaled_fourier(x, alpha=1):
    x = np.asarray(x)
    return czt(x, np.exp(2j*np.pi*alpha/x.size))


def scaled_fourier2d(x, alpha=1):
    x = np.asarray(x)
    if x.ndim != 2:
        raise ValueError('Use a two-dimensional array.')
    y = np.apply_along_axis(scaled_fourier, 0, x, alpha)
    return np.apply_along_axis(scaled_fourier, 1, y, alpha)


def toeplitz_product(column, row, x):
    column, row, x = map(np.asarray, (column, row, x))
    m, n = column.size, row.size
    if min(m, n) == 0 or row[0] != column[0] or x.shape != (n,):
        raise ValueError('Use compatible first row, first column, and vector.')
    size = 1 << (m+n-2).bit_length()
    kernel = np.concatenate((column, np.zeros(size-m-n+1), row[:0:-1]))
    y = np.fft.ifft(np.fft.fft(kernel)*np.fft.fft(x, size))[:m]
    return y.real if all(np.isrealobj(v) for v in (column, row, x)) else y


def inverse_modulo(x, modulus):
    return pow(int(x), -1, int(modulus))


def rader_fft(x, generator):
    """Prime-length DFT via a cyclic convolution of length p-1."""
    from sympy import isprime
    x = np.asarray(x, dtype=complex)
    p = x.size
    indices = np.array([pow(int(generator), k, p) for k in range(p-1)])
    if not isprime(p) or set(indices) != set(range(1, p)):
        raise ValueError('Use prime length and a generator of its nonzero residues.')
    inverse = np.array([inverse_modulo(j, p) for j in indices])
    kernel = np.exp(-2j*np.pi*inverse/p)
    convolution = np.fft.ifft(np.fft.fft(x[indices])*np.fft.fft(kernel))
    y = np.empty(p, dtype=complex)
    y[0], y[inverse] = x.sum(), x[0]+convolution
    return y


def simpson_fourier(samples, start, step, frequency, spacing, count):
    """Approximate integral f(t) exp(-i*xi*t) dt at xi=frequency+k*spacing."""
    samples = np.asarray(samples)
    n = samples.size
    if n < 3 or n % 2 != 1 or step <= 0:
        raise ValueError('Simpson quadrature needs an odd number of at least 3 samples.')
    weights = np.ones(n)
    weights[1:-1:2], weights[2:-1:2] = 4, 2
    j = np.arange(n)
    x = weights*samples*np.exp(-1j*step*frequency*j)
    xi = frequency+spacing*np.arange(count)
    values = czt(x, np.exp(1j*step*spacing), count)
    return (step/3)*np.exp(-1j*start*xi)*values


def zoom_spectrum(x, start, spacing, count):
    """Finite-record spectrum at start+k*spacing cycles per original record."""
    x = np.asarray(x)
    modulated = x*np.exp(-2j*np.pi*start*np.arange(x.size)/x.size)
    return czt(modulated, np.exp(2j*np.pi*spacing/x.size), count)
