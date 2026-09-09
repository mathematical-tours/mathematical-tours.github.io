"""Educational DFT algorithms using the negative-exponential convention."""
import numpy as np
from math import gcd


def phase_multiply(x, alpha):
    x = np.asarray(x, dtype=complex)
    if x.ndim != 1 or x.size == 0:
        raise ValueError('Expected a nonempty vector.')
    return x * np.exp(-2j*np.pi*alpha*np.arange(x.size)/x.size)


def fft_recursive(x, inverse=False):
    x = np.asarray(x, dtype=complex)
    n = x.size
    if x.ndim != 1 or n == 0 or n & (n-1):
        raise ValueError('Expected a vector of power-of-two length.')
    if inverse:
        return np.conj(fft_recursive(np.conj(x)))/n
    if n == 1:
        return x.copy()
    a = fft_recursive(x[::2])
    b = phase_multiply(fft_recursive(x[1::2]), .5)
    return np.r_[a+b, a-b]


def bit_reversal_indices(n):
    if n < 1 or n & (n-1):
        raise ValueError('Expected a positive power of two.')
    indices = np.array([0])
    while indices.size < n:
        indices = np.r_[2*indices, 2*indices+1]
    return indices


def reverse_index(i, bits):
    if bits < 0 or not 0 <= i < 2**bits:
        raise ValueError('Index does not fit in the specified bit count.')
    result = 0
    for _ in range(bits):
        result, i = 2*result+(i & 1), i >> 1
    return result


def bit_reverse_recursive(x):
    x = np.asarray(x)
    n = x.size
    if x.ndim != 1 or n == 0 or n & (n-1):
        raise ValueError('Expected a vector of power-of-two length.')
    if n == 1:
        return x.copy()
    return np.r_[bit_reverse_recursive(x[::2]), bit_reverse_recursive(x[1::2])]


def fft_dit(x, inverse=False):
    x = np.asarray(x, dtype=complex)
    if x.ndim != 1:
        raise ValueError('Expected a vector.')
    n = x.size
    y = x[bit_reversal_indices(n)].copy()
    size = 2
    while size <= n:
        half = size//2
        w = np.exp((-1 if not inverse else 1)*2j*np.pi*np.arange(half)/size)
        for start in range(0, n, size):
            a = y[start:start+half].copy()
            b = w*y[start+half:start+size]
            y[start:start+half], y[start+half:start+size] = a+b, a-b
        size *= 2
    return y/n if inverse else y


def fft_dif(x, inverse=False):
    y = np.array(x, dtype=complex, copy=True)
    if y.ndim != 1:
        raise ValueError('Expected a vector.')
    n, size = y.size, y.size
    indices = bit_reversal_indices(n)
    while size > 1:
        half = size//2
        w = np.exp((-1 if not inverse else 1)*2j*np.pi*np.arange(half)/size)
        for start in range(0, n, size):
            a, b = y[start:start+half].copy(), y[start+half:start+size].copy()
            y[start:start+half], y[start+half:start+size] = a+b, w*(a-b)
        size //= 2
    return y[indices]/n if inverse else y[indices]


def good_thomas(x, p):
    x = np.asarray(x, dtype=complex)
    if x.ndim != 1 or p < 2:
        raise ValueError('Expected a vector and a factor p >= 2.')
    n, q = x.size, x.size//p
    if p < 2 or q < 2 or p*q != n or gcd(p, q) != 1:
        raise ValueError('Expected coprime factors p and q, both at least two.')
    F = np.empty((p, q), dtype=complex)
    indices = np.arange(n)
    F[indices % p, indices % q] = x
    spectrum = np.fft.fft2(F)
    s1, s2 = np.indices((p, q))
    y = np.empty(n, dtype=complex)
    y[(s1*q+s2*p) % n] = spectrum
    return y


def cooley_tukey(x, p):
    x = np.asarray(x, dtype=complex)
    if x.ndim != 1:
        raise ValueError('Expected a vector.')
    n = x.size
    if p < 2 or n % p or n//p < 2:
        raise ValueError('Expected a nontrivial factor p of the vector length.')
    q = n//p
    F = np.fft.fft(x.reshape(q, p), axis=0)
    d, b = np.indices((q, p))
    F *= np.exp(-2j*np.pi*d*b/n)
    return np.fft.fft(F, axis=1).T.ravel()


def split_radix(x):
    x = np.asarray(x, dtype=complex)
    n = x.size
    if x.ndim != 1 or n == 0 or n & (n-1):
        raise ValueError('Expected a vector of power-of-two length.')
    if n == 1:
        return x.copy()
    if n == 2:
        return np.array([x[0]+x[1], x[0]-x[1]])
    E = split_radix(x[::2])
    w = np.exp(-2j*np.pi*np.arange(n//4)/n)
    u, v = w*split_radix(x[1::4]), w**3*split_radix(x[3::4])
    a, b = u+v, -1j*(u-v)
    return np.r_[E[:n//4]+a, E[n//4:]+b, E[:n//4]-a, E[n//4:]-b]


def overlap_add(f, g):
    f, g = np.asarray(f), np.asarray(g)
    N, M = f.size, g.size
    if f.ndim != 1 or g.ndim != 1 or min(N, M) == 0:
        raise ValueError('Expected two nonempty vectors.')
    L = 1 << (2*M-2).bit_length()
    spectrum = np.fft.fft(g, L)
    result = np.zeros(N+M-1, dtype=complex)
    for start in range(0, N, M):
        block = f[start:start+M]
        values = np.fft.ifft(np.fft.fft(block, L)*spectrum)
        length = block.size+M-1
        result[start:start+length] += values[:length]
    return np.real_if_close(result) if np.isrealobj(f) and np.isrealobj(g) else result
