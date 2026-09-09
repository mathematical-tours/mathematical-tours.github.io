"""Polynomial products and base-b integer arithmetic, with exact checks."""
import numpy as np


def polynomial_add(p, q):
    n = max(len(p), len(q))
    return [(p[k] if k < len(p) else 0)+(q[k] if k < len(q) else 0) for k in range(n)]


def karatsuba(p, q):
    size = len(p)+len(q)-1
    if not p or not q:
        return []
    n = 1 << (max(len(p), len(q))-1).bit_length()
    p, q = list(p)+[0]*(n-len(p)), list(q)+[0]*(n-len(q))
    if n == 1:
        return [p[0]*q[0]]
    k = n//2
    r0, r2 = karatsuba(p[:k], q[:k]), karatsuba(p[k:], q[k:])
    r1 = karatsuba(polynomial_add(p[:k], p[k:]), polynomial_add(q[:k], q[k:]))
    r1 = polynomial_add(polynomial_add(r1, [-v for v in r0]), [-v for v in r2])
    result = polynomial_add(polynomial_add(r0, [0]*k+r1), [0]*(2*k)+r2)
    return result[:size]


def integer_to_digits(value, base):
    if not isinstance(value, int) or value < 0 or not isinstance(base, int) or base < 2:
        raise ValueError('Use a nonnegative integer and an integer base >= 2.')
    digits = []
    while value:
        value, digit = divmod(value, base)
        digits.append(digit)
    return digits or [0]


def digits_to_integer(digits, base):
    if not isinstance(base, int) or base < 2:
        raise ValueError('Base must be an integer >= 2.')
    value = 0
    for digit in reversed(digits):
        if int(digit) != digit or not 0 <= digit < base:
            raise ValueError('Digits must be integers in range(base).')
        value = base*value+int(digit)
    return value


def multiply_digits(x, y, base):
    exact = digits_to_integer(x, base)*digits_to_integer(y, base)
    if not x or not y or min(len(x), len(y))*(base-1)**2 > 2**40:
        raise ValueError('Use nonempty digit vectors with a moderate coefficient bound.')
    length = 1 << (len(x)+len(y)-2).bit_length()
    product = np.fft.irfft(np.fft.rfft(x, length)*np.fft.rfft(y, length), length)
    coefficients = [int(v) for v in np.rint(product[:len(x)+len(y)-1])]
    if min(coefficients) < 0:
        raise ArithmeticError('FFT error produced a negative coefficient.')
    digits, carry = [], 0
    for coefficient in coefficients:
        carry, digit = divmod(coefficient+carry, base)
        digits.append(digit)
    while carry:
        carry, digit = divmod(carry, base)
        digits.append(digit)
    while len(digits) > 1 and digits[-1] == 0:
        digits.pop()
    if digits_to_integer(digits, base) != exact:
        raise ArithmeticError('FFT product failed exact verification.')
    return digits


def integer_product_example():
    a, b, base = 262122154512154212, 314134464653513212, 20
    digits = multiply_digits(integer_to_digits(a, base), integer_to_digits(b, base), base)
    return digits_to_integer(digits, base) == a*b
