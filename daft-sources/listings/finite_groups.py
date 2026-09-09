"""Finite-group examples. NumPy uses the negative-exponential DFT."""
import numpy as np
import sympy as sp


def circulant_determinant(first_row):
    eigenvalues = np.fft.fft(np.asarray(first_row, dtype=complex))
    return np.prod(eigenvalues)


def rotation_product_trace():
    a = sp.symbols('alpha', real=True)
    c, s = sp.cos(a), sp.sin(a)
    rx = sp.Matrix([[1, 0, 0], [0, c, -s], [0, s, c]])
    ry = sp.Matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return sp.simplify(sp.trace(rx * ry.T))


def random_walk(p0, v, k):
    """Law after k >= 0 steps; p0 and v are probability vectors of equal size."""
    spectrum = np.fft.fft(p0) * np.fft.fft(v)**k
    return np.fft.ifft(spectrum).real


def paley_hadamard(p):
    if p % 4 != 3 or not sp.isprime(p):
        raise ValueError('Expected a prime congruent to 3 modulo 4.')
    Q = np.array([[int(sp.legendre_symbol((j-i) % p, p))
                   for j in range(p)] for i in range(p)])
    H = np.ones((p+1, p+1), dtype=int)
    H[1:, 1:] = Q - np.eye(p, dtype=int)
    return H


def paley_example():
    H = paley_hadamard(31)
    return H @ H.T
