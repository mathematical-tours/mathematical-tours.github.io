"""Small exact coding examples; vector coordinates increase with polynomial degree."""
import numpy as np
import sympy as sp
import galois
from transforms import fwt
from finite_fields import field_with_root, finite_dft


def binary_digits(x, width):
    if not 0 <= x < 2**width:
        raise ValueError('The nonnegative integer must fit the requested width.')
    return np.array([(x >> j) & 1 for j in range(width)], dtype=int)


def encode_rm(a, b, variables):
    if not 0 <= a < 2**variables or b not in (0, 1):
        raise ValueError('Use a valid linear part and binary constant.')
    impulse = np.zeros(2**variables, dtype=int)
    impulse[a] = 1
    signs = fwt(impulse).real.astype(int)
    return ((1-signs)//2 + b) % 2


def decode_rm(word):
    word = np.asarray(word)
    n = word.size
    if word.ndim != 1 or n < 2 or n & (n-1) or not np.all((word == 0) | (word == 1)):
        raise ValueError('Use a binary vector of dyadic length at least two.')
    spectrum = fwt(1-2*word.astype(int)).real
    a = int(np.argmax(abs(spectrum)))
    b = int(spectrum[a] < 0)
    return encode_rm(a, b, n.bit_length()-1)


def hamming_check(bits):
    if bits < 2:
        raise ValueError('Use at least two parity-check coordinates.')
    return np.array([binary_digits(j, bits) for j in range(1, 2**bits)]).T


def decode_hamming(word):
    word = np.asarray(word)
    n = word.size
    if word.ndim != 1 or n < 3 or (n+1) & n or not np.all((word == 0) | (word == 1)):
        raise ValueError('Use a binary Hamming word of length 2^k-1, k>=2.')
    bits = (n+1).bit_length()-1
    syndrome = hamming_check(bits) @ word.astype(int) % 2
    position = sum(int(bit) << j for j, bit in enumerate(syndrome))
    corrected = word.astype(int).copy()
    if position:
        corrected[position-1] ^= 1
    return corrected


def hamming_weight_distribution(bits=3):
    x = sp.Symbol('Y')
    n = 2**bits-1
    polynomial = sp.Poly(sp.expand(((1+x)**n+n*(1-x)**(2**(bits-1))
                                   *(1+x)**(2**(bits-1)-1))/2**bits), x)
    return [int(polynomial.nth(k)) for k in range(n+1)]


def cyclotomic_coset(p, n, index):
    orbit = []
    while index not in orbit:
        orbit.append(index)
        index = p*index % n
    return orbit


def bch_code(p=2, n=31, distance=7, modulus=None):
    if not 2 <= distance <= n:
        raise ValueError('The designed distance must lie between 2 and n.')
    field, root = field_with_root(p, n, modulus)
    indices = sorted({j for k in range(1, distance) for j in cyclotomic_coset(p, n, k)})
    polynomial = galois.Poly.Roots(root**np.array(indices))
    assert np.all(polynomial.coeffs**p == polynomial.coeffs)
    prime = galois.GF(p, compile='python-calculate')
    generator = prime([int(c) for c in polynomial.coeffs[::-1]])
    return dict(p=p, n=n, k=n-polynomial.degree, t=(distance-1)//2,
                field=field, root=root, generator=generator, distance=distance)


def encode_bch(message, code):
    prime = type(code['generator'])
    message = prime(message)
    if message.shape != (code['k'],):
        raise ValueError('The message length must equal the code dimension.')
    word = prime.Zeros(code['n'])
    g = code['generator']
    for j, value in enumerate(message):
        word[j:j+g.size] += value*g
    return word


def error_locator(syndromes):
    """Minimal locator from S_1,...,S_(2t), including fewer-than-t errors."""
    field = type(syndromes)
    t = syndromes.size//2
    indices = np.add.outer(np.arange(t), np.arange(t))
    hankel = syndromes[indices]
    degree = int(np.count_nonzero(np.any(hankel.row_reduce() != 0, axis=1)))
    if degree == 0:
        if np.any(syndromes):
            raise ValueError('The syndromes are inconsistent with the error bound.')
        return field([1])
    indices = degree+np.arange(degree)[:, None]-np.arange(1, degree+1)[None, :]
    tail = np.linalg.solve(syndromes[indices], -syndromes[degree:2*degree])
    locator = field(np.r_[1, tail])
    for j in range(degree, 2*t):
        if syndromes[j]+np.sum(tail*syndromes[j-np.arange(1, degree+1)]) != 0:
            raise ValueError('No locator fits the supplied syndromes.')
    return locator


def decode_bch(received, code, method='spectrum'):
    field, root, n, t = code['field'], code['root'], code['n'], code['t']
    prime = type(code['generator'])
    received = prime(received)
    if received.shape != (n,) or t < 1:
        raise ValueError('Use a code with positive correction radius and matching word length.')
    S = (root**np.outer(np.arange(1, 2*t+1), np.arange(n))) @ field(received)
    locator = error_locator(S)
    degree = locator.size-1
    if method == 'spectrum':
        spectrum = field.Zeros(n+1)
        spectrum[1:2*t+1] = S
        for j in range(2*t+1, n+1):
            if degree:
                spectrum[j] = -np.sum(locator[1:]*spectrum[j-np.arange(1, degree+1)])
        spectrum[0] = spectrum[n]
        errors = finite_dft(spectrum[:n], root**-1, inverse=True)
    elif method == 'roots':
        values = (root**(-np.outer(np.arange(n), np.arange(degree+1)))) @ locator
        positions = np.flatnonzero(values == 0)
        if positions.size != degree:
            raise ValueError('The locator does not have the required number of roots.')
        errors = field.Zeros(n)
        if degree:
            matrix = root**np.outer(np.arange(1, degree+1), positions)
            errors[positions] = np.linalg.solve(matrix, S[:degree])
    else:
        raise ValueError('Choose spectrum or roots.')
    if np.any(errors**code['p'] != errors) or np.count_nonzero(errors) > t:
        raise ValueError('The reconstructed error is outside the permitted alphabet or radius.')
    corrected = field(received)-errors
    check = (root**np.outer(np.arange(1, code['distance']), np.arange(n))) @ corrected
    if np.any(check):
        raise ValueError('The corrected word does not satisfy the BCH checks.')
    return prime([int(c) for c in corrected])


def bch_example():
    code = bch_code()
    message = np.random.default_rng(61).integers(0, 2, code['k'])
    word = encode_bch(message, code)
    received = word.copy()
    received[[0, 12, 30]] += type(word)(1)
    corrected = decode_bch(received, code)
    assert np.array_equal(corrected, word)
    assert np.array_equal(corrected, decode_bch(received, code, method='roots'))
    return code, word, received, corrected
