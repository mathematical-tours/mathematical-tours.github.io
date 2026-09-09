"""Elementary fast transforms; indices start at zero.

FWT and tensor transforms are unnormalized. Haar analysis is orthonormal and
packs fine-to-coarse details followed by the final scaling coefficient.
"""
import numpy as np


def fwt(x):
    x = np.asarray(x, dtype=complex)
    n = x.size
    if x.ndim != 1 or n == 0 or n & (n-1):
        raise ValueError('Expected a nonempty vector of power-of-two length.')
    if n == 1:
        return x.copy()
    u, v = fwt(x[:n//2]), fwt(x[n//2:])
    return np.r_[u + v, u - v]


def fwt2d(x):
    rows = np.apply_along_axis(fwt, 1, x)
    return np.apply_along_axis(fwt, 0, rows)


def haar(x):
    x = np.asarray(x, dtype=complex)
    n = x.size
    if x.ndim != 1 or n == 0 or n & (n-1):
        raise ValueError('Expected a nonempty vector of power-of-two length.')
    details = []
    while x.size > 1:
        even, odd = x[::2], x[1::2]
        details.append((even - odd) / np.sqrt(2))
        x = (even + odd) / np.sqrt(2)
    return np.concatenate(details + [x])


def haar_inverse(y):
    y = np.asarray(y, dtype=complex)
    n = y.size
    if y.ndim != 1 or n == 0 or n & (n-1):
        raise ValueError('Expected a nonempty vector of power-of-two length.')
    x, size = y[-1:].copy(), 1
    while size < n:
        detail = y[n-2*size:n-size]
        z = np.empty(2*size, dtype=complex)
        z[::2], z[1::2] = (x+detail)/np.sqrt(2), (x-detail)/np.sqrt(2)
        x, size = z, 2*size
    return x


def walsh_sequency(n):
    if n < 1 or n & (n-1):
        raise ValueError('Expected a positive power of two.')
    counts, parity = np.array([0]), np.array([0])
    for _ in range(n.bit_length()-1):
        parity = np.r_[parity, 1-parity]
        counts = 2*np.tile(counts, 2) + parity
    return counts


def tensor_transform(x, A):
    x, A = np.asarray(x, dtype=complex), np.asarray(A, dtype=complex)
    s = A.shape[0]
    if A.shape != (s, s) or s < 2 or x.ndim != 1 or x.size == 0:
        raise ValueError('Expected a square base matrix and a nonempty vector.')
    if x.size == 1:
        return x.copy()
    if x.size % s:
        raise ValueError('Vector length must be a power of the base size.')
    blocks = np.stack([tensor_transform(part, A) for part in np.split(x, s)])
    return (A @ blocks).ravel()


def walsh_approx2d(x, M):
    x = np.asarray(x)
    if x.ndim != 2 or not 0 <= M <= x.size:
        raise ValueError('Expected an image and a valid coefficient count.')
    sy, sx = walsh_sequency(x.shape[0]), walsh_sequency(x.shape[1])
    vertical, horizontal = np.meshgrid(sy, sx, indexing='ij')
    order = np.lexsort((vertical.ravel(), (vertical+horizontal).ravel()))
    spectrum = fwt2d(x).ravel()
    spectrum[order[M:]] = 0
    result = fwt2d(spectrum.reshape(x.shape)) / x.size
    return np.real_if_close(result) if np.isrealobj(x) else result
