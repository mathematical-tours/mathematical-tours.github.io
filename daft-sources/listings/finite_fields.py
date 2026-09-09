"""Exact Fourier transforms over finite fields and integer residue rings."""
from math import gcd
import numpy as np
import sympy as sp
import galois
from polynomials import integer_to_digits, digits_to_integer


def cyclotomic_example(n=105):
    x = sp.Symbol('X')
    return sp.cyclotomic_poly(n, x, polys=True)


def cyclotomic_coefficients(n):
    """Ascending exact coefficients from the truncated Mobius product."""
    if n < 1:
        raise ValueError('The cyclotomic index must be positive.')
    if n == 1:
        return [-1, 1]
    degree = int(sp.totient(n))
    coefficients = [1]+[0]*degree
    for divisor in sp.divisors(n):
        sign = sp.mobius(n//divisor)
        if sign == 1:
            for k in range(degree, divisor-1, -1):
                coefficients[k] -= coefficients[k-divisor]
        elif sign == -1:
            for k in range(divisor, degree+1):
                coefficients[k] += coefficients[k-divisor]
    return coefficients


def first_cyclotomic_coefficients(largest=7, limit=3200):
    """First positive index with a coefficient of each absolute value 1,...,largest."""
    if largest < 1 or limit < 1:
        raise ValueError('Use positive bounds.')
    found = {1: 1}
    if largest == 1:
        return found
    for n in range(3, limit+1, 2):
        factors = sp.factorint(n)
        if any(exponent != 1 for exponent in factors.values()):
            continue  # Repeated prime factors only insert zeros; 2 changes signs.
        if len(factors) < 3:
            continue  # One or two distinct primes give only 0, 1, and -1.
        coefficients = cyclotomic_coefficients(n)
        for value in set(abs(int(c)) for c in coefficients):
            if 1 <= value <= largest:
                found.setdefault(value, n)
        if len(found) == largest:
            return found
    return found  # Missing values have not been found within the explicit limit.


def field_with_root(p, n, modulus=None):
    """A minimal extension of F_p containing a root of exact order n."""
    if not sp.isprime(p) or n < 1 or gcd(p, n) != 1:
        raise ValueError('Use prime p and positive n coprime to p.')
    prime = galois.GF(p, compile='python-calculate')
    if modulus is None:
        coefficients = [int(c) % p for c in cyclotomic_example(n).all_coeffs()]
        factors, _ = galois.Poly(coefficients, field=prime).factors()
        polynomial = factors[0]
    else:
        polynomial = galois.Poly(modulus, field=prime)
    if polynomial.degree == 1:
        field, root = prime, -polynomial.coeffs[1]/polynomial.coeffs[0]
    else:
        field = galois.GF(p**polynomial.degree, irreducible_poly=polynomial,
                          compile='python-calculate')
        root = field(p)  # The class of X in the polynomial-basis representation.
    if root**n != 1 or any(root**(n//ell) == 1 for ell in sp.factorint(n)):
        raise ValueError('The selected polynomial does not give a root of order n.')
    return field, root


def finite_dft(x, root, inverse=False):
    field = type(root)
    x = field(x)
    n = x.size
    if x.ndim != 1 or n == 0 or root**n != 1:
        raise ValueError('The input length must match the root order.')
    if any(root**(n//ell) == 1 for ell in sp.factorint(n)):
        raise ValueError('The root must have exact order n.')
    powers = np.outer(np.arange(n), np.arange(n))
    y = (root**(powers if inverse else -powers)) @ x
    return y/field(n % field.characteristic) if inverse else y


def finite_fft(x, root):
    """Mixed-radix finite-field FFT, with a direct sum at prime leaves."""
    field = type(root)
    x = field(x)
    n = x.size
    if n == 1:
        return x.copy()
    radix = min(sp.factorint(n))
    if radix == n:
        return finite_dft(x, root)
    m = n//radix
    blocks = [finite_fft(x[j::radix], root**radix) for j in range(radix)]
    k = np.arange(n)
    y = field.Zeros(n)
    for j in range(radix):
        y += root**(-j*k)*blocks[j][k % m]
    return y


def finite_transform_example():
    field, root = field_with_root(2, 15)
    x = field(np.arange(15) % 2)
    y = finite_dft(x, root)
    assert np.array_equal(finite_dft(y, root, inverse=True), x)
    assert np.array_equal(finite_fft(x, root), y)
    return x, y


def _modular_fft(x, root, modulus):
    n = len(x)
    if n == 1:
        return x.copy()
    even = _modular_fft(x[::2], root*root % modulus, modulus)
    odd = _modular_fft(x[1::2], root*root % modulus, modulus)
    inverse_root = pow(root, -1, modulus)
    mixed = [value*pow(inverse_root, k, modulus) % modulus for k, value in enumerate(odd)]
    return ([(a+b) % modulus for a, b in zip(even, mixed)]
            + [(a-b) % modulus for a, b in zip(even, mixed)])


def modular_fft(x, modulus, root, inverse=False):
    x = [int(a) % modulus for a in x]
    n = len(x)
    if modulus < 2 or n == 0 or n & (n-1) or gcd(n, modulus) != 1:
        raise ValueError('Use dyadic length invertible modulo the modulus.')
    if pow(root, n, modulus) != 1 or (n > 1 and pow(root, n//2, modulus) != modulus-1):
        raise ValueError('The root must be a principal dyadic root.')
    y = _modular_fft(x, pow(root, -1, modulus) if inverse else root, modulus)
    return [a*pow(n, -1, modulus) % modulus for a in y] if inverse else y


def multiply_integers_modular(a, b, base=10000):
    """Exact integer multiplication via an NTT, without floating-point rounding."""
    if a == 0 or b == 0:
        return 0
    x, y = integer_to_digits(abs(a), base), integer_to_digits(abs(b), base)
    size = max(2, 1 << (len(x)+len(y)-2).bit_length())
    bound = min(len(x), len(y))*(base-1)**2
    step = max(1, (bound.bit_length()+size//2-1)//(size//2))
    modulus, root = (1 << (step*size//2))+1, 1 << step
    X = modular_fft(x+[0]*(size-len(x)), modulus, root)
    Y = modular_fft(y+[0]*(size-len(y)), modulus, root)
    coefficients = modular_fft([u*v % modulus for u, v in zip(X, Y)], modulus, root, True)
    digits, carry = [], 0
    for value in coefficients:
        carry, digit = divmod(value+carry, base)
        digits.append(digit)
    while carry:
        carry, digit = divmod(carry, base)
        digits.append(digit)
    return (1 if (a > 0) == (b > 0) else -1)*digits_to_integer(digits, base)
