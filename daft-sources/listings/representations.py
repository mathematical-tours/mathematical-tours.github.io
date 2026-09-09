"""Affine-group wavelets, with a negative-phase ordinary forward DFT."""
import numpy as np
import sympy as sp


def finite_wavelet(f, psi):
    f, psi = np.asarray(f, dtype=complex), np.asarray(psi, dtype=complex)
    p = psi.size
    if psi.ndim != 1 or f.shape != (p,) or not sp.isprime(p):
        raise ValueError('Use two vectors of the same prime length.')
    F, P = np.fft.fft(f), np.fft.fft(psi)
    k = np.arange(p)
    return np.array([np.fft.ifft(F*np.conj(P[a*k % p])) for a in range(1, p)])


def inverse_finite_wavelet(coefficients, psi, mode='mean-zero'):
    psi = np.asarray(psi, dtype=complex)
    p = psi.size
    coefficients = np.asarray(coefficients, dtype=complex)
    if psi.ndim != 1 or not sp.isprime(p) or coefficients.shape != (p-1, p):
        raise ValueError('Use p-1 dilation rows and p translation columns.')
    P, k = np.fft.fft(psi), np.arange(p)
    c = np.sum(abs(P[1:])**2)
    if c <= 0:
        raise ValueError('The wavelet needs a nonzero nonconstant component.')
    if mode == 'mean-zero':
        if abs(P[0]) > 1e-10*np.linalg.norm(P):
            raise ValueError('The mean-zero mode requires a mean-zero wavelet.')
    elif mode == 'tight':
        if not np.isclose((p-1)*abs(P[0])**2, c, rtol=1e-10, atol=0):
            raise ValueError('The full-space mode requires the stated tight-frame condition.')
    else:
        raise ValueError('Choose mean-zero or tight.')
    spectrum = np.zeros(p, dtype=complex)
    for a in range(1, p):
        spectrum += np.fft.fft(coefficients[a-1])*P[a*k % p]
    if mode == 'mean-zero':
        spectrum[0] = 0  # The analysis contains no information about the signal mean.
    return np.fft.ifft(spectrum/c)
