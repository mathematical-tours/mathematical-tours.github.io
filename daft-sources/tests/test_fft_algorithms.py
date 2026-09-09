"""Compare educational FFTs with explicit Fourier sums and direct convolution."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'listings'))
from fft_algorithms import (fft_recursive, fft_dit, fft_dif, split_radix,
    bit_reversal_indices, reverse_index, bit_reverse_recursive,
    phase_multiply, good_thomas, cooley_tukey, overlap_add)


class FourierAlgorithms(unittest.TestCase):
    def test_radix_two_and_split_radix(self):
        rng = np.random.default_rng(100)
        for n in [1, 2, 4, 8, 16, 32, 64]:
            x = rng.normal(size=n)+1j*rng.normal(size=n)
            matrix = np.exp(-2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)
            expected = matrix @ x
            for method in [fft_recursive, fft_dit, fft_dif, split_radix]:
                np.testing.assert_allclose(method(x), expected, atol=3e-13)
            for method in [fft_recursive, fft_dit, fft_dif]:
                np.testing.assert_allclose(method(x,inverse=True), matrix.conj()@x/n, atol=1e-14)
            saved = x.copy()
            for method in [fft_recursive, fft_dit, fft_dif, split_radix]:
                method(x)
                np.testing.assert_array_equal(x, saved)

    def test_bit_reversal(self):
        for bits in range(9):
            n = 2**bits
            expected = [int(format(i, f'0{bits}b')[::-1], 2) for i in range(n)]
            np.testing.assert_array_equal(bit_reversal_indices(n), expected)
            np.testing.assert_array_equal(bit_reverse_recursive(np.arange(n)), expected)
            for i in range(n): self.assertEqual(reverse_index(i,bits), expected[i])

    def test_phase_operator(self):
        x = np.arange(7)+1j
        for alpha in [-1.3, 0, .5, 3]:
            expected = np.array([x[j]*np.exp(-2j*np.pi*alpha*j/7) for j in range(7)])
            np.testing.assert_allclose(phase_multiply(x,alpha), expected)

    def test_mixed_radix_and_good_thomas(self):
        rng = np.random.default_rng(71)
        for p,q in [(2,3),(3,4),(4,5),(5,7),(6,5),(3,3),(4,4),(2,8)]:
            n = p*q;x = rng.normal(size=n)+1j*rng.normal(size=n)
            expected = np.exp(-2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)@x
            np.testing.assert_allclose(cooley_tukey(x,p),expected,atol=2e-13)
            if np.gcd(p,q)==1:
                np.testing.assert_allclose(good_thomas(x,p),expected,atol=2e-13)

    def test_overlap_add_including_tail_and_scalar_kernel(self):
        rng = np.random.default_rng(73)
        for n,m in [(1,1),(2,5),(21,3),(23,5),(67,9),(13,1)]:
            f = rng.normal(size=n)+1j*rng.normal(size=n)
            g = rng.normal(size=m)+1j*rng.normal(size=m)
            np.testing.assert_allclose(overlap_add(f,g),np.convolve(f,g),atol=4e-15)
            np.testing.assert_allclose(overlap_add(f.real,g.real),np.convolve(f.real,g.real),atol=4e-15)


if __name__ == '__main__': unittest.main()
