"""Interpolation, differentiation, and the spectral structure of a DFT."""
import numpy as np


def trigonometric_interpolate(x, eta):
    x = np.asarray(x)
    n = x.size
    if x.ndim != 1 or n % 2 != 1 or not isinstance(eta, int) or eta < 1:
        raise ValueError('Expected an odd-length vector and an integer factor.')
    spectrum, middle = np.fft.fft(x), (n+1)//2
    padded = eta*np.r_[spectrum[:middle], np.zeros((eta-1)*n), spectrum[middle:]]
    y = np.fft.ifft(padded)
    return np.real_if_close(y) if np.isrealobj(x) else y


def dct2(x):
    x = np.asarray(x)
    if x.ndim != 1 or x.size == 0 or not np.isrealobj(x):
        raise ValueError('Expected a nonempty real vector.')
    n = x.size
    padded = np.zeros(4*n)
    padded[1:2*n:2] = x
    return np.fft.fft(padded).real[:n]


def dct3(x):
    x = np.asarray(x)
    if x.ndim != 1 or x.size == 0 or not np.isrealobj(x):
        raise ValueError('Expected a nonempty real vector.')
    n = x.size
    return np.fft.fft(x, 4*n).real[1:2*n:2]-x[0]/2


def chebyshev_interpolate(values, points):
    coefficients = 2*dct2(values)/len(values)
    coefficients[0] /= 2
    return np.polynomial.chebyshev.chebval(points, coefficients)


def chebyshev_example(n=16, alpha=.3):
    nodes = np.cos((np.arange(n)+.5)*np.pi/n)
    values = 1/(alpha**2+nodes**2)
    grid = np.linspace(-1, 1, 201)
    return nodes, values, grid, chebyshev_interpolate(values, grid)


def fractional_derivative(f, alpha, period):
    f = np.asarray(f)
    n = f.size
    if f.ndim != 1 or n % 2 != 1 or alpha < 0 or period <= 0:
        raise ValueError('Use odd sample count, nonnegative order, positive period.')
    omega = 2*np.pi*np.fft.fftfreq(n, d=period/n)
    multiplier = np.ones(n, dtype=complex)
    if alpha > 0:
        multiplier = abs(omega)**alpha*np.exp(.5j*np.pi*alpha*np.sign(omega))
    y = np.fft.ifft(multiplier*np.fft.fft(f))
    return y.real if np.isrealobj(f) else y


def fourier_components(x):
    x = np.asarray(x, dtype=complex)
    n = x.size
    if x.ndim != 1 or n == 0:
        raise ValueError('Expected a nonempty vector.')
    transformed = np.fft.fft(x)/np.sqrt(n)
    reverse = (-np.arange(n)) % n
    powers = np.stack([x, transformed, x[reverse], transformed[reverse]])
    eigenvalues = np.array([1, -1, 1j, -1j])
    weights = eigenvalues[:, None]**(-np.arange(4)[None, :])/4
    return weights @ powers


def fractional_dft(x, alpha, turns=(0, 0, 0, 0)):
    turns = np.asarray(turns)
    if turns.shape != (4,) or np.any(turns != np.round(turns)):
        raise ValueError('Choose one integer branch shift per eigenspace.')
    phases = np.array([0, np.pi, np.pi/2, -np.pi/2])+2*np.pi*turns
    return np.exp(1j*alpha*phases) @ fourier_components(x)


def fourier_eigenbasis(n):
    if n < 3:
        raise ValueError('This construction uses n >= 3.')
    eye = np.eye(n)
    S = np.roll(eye, 1, axis=0)+np.roll(eye, -1, axis=0)
    S += np.diag(2*np.cos(2*np.pi*np.arange(n)/n)-4)
    even, odd = [eye[:, 0]], []
    for k in range(1, (n+1)//2):
        even.append((eye[:, k]+eye[:, -k])/np.sqrt(2))
        odd.append((eye[:, k]-eye[:, -k])/np.sqrt(2))
    if n % 2 == 0:
        even.append(eye[:, n//2])
    vectors, values, parity = [], [], []
    for sign, basis in enumerate([even, odd]):
        B = np.column_stack(basis)
        eigenvalues, V = np.linalg.eigh(B.T @ S @ B)
        vectors.append(B @ V); values.extend(eigenvalues)
        parity.extend([sign]*len(eigenvalues))
    order = np.lexsort((parity, -np.round(values, 12)))
    V = np.column_stack(vectors)[:, order]
    pivots = np.argmax(abs(V), axis=0)
    V *= np.sign(V[pivots, np.arange(n)])
    return V


def orthogonalize_translates(x):
    x = np.asarray(x)
    spectrum = np.fft.fft(x)
    if np.any(abs(spectrum) == 0):
        raise ValueError('Every Fourier coefficient must be nonzero.')
    y = np.fft.ifft(np.sqrt(x.size)*spectrum/abs(spectrum))
    return np.real_if_close(y) if np.isrealobj(x) else y
