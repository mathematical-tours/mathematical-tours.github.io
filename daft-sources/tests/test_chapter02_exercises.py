"""Independent checks for the repaired duality exercises."""
import os
import sys
import unittest
from itertools import product
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('NUMBA_CACHE_DIR', str(ROOT / 'build/.numba-cache'))
sys.path.insert(0, str(ROOT / 'listings'))
from transforms import fwt, walsh_approx2d


class DualityExercises(unittest.TestCase):
    def test_reciprocity_lattice_pairing(self):
        for p, r in product([3, 5, 7, 11, 13, 17, 19, 23], repeat=2):
            if p == r:
                continue
            even = sum(r*a//p for a in range(2, p, 2))
            below = sum(r*a//p for a in range(1, (p+1)//2))
            above = sum(p*b//r for b in range(1, (r+1)//2))
            self.assertEqual(even % 2, below % 2)
            self.assertEqual(below+above, (p-1)*(r-1)//4)
            self.assertEqual((-1)**below, int(sp.legendre_symbol(r, p)))

    def test_factorial_effects(self):
        x = np.array([69, 81, 63, 77, 61, 92, 54, 89])
        expected = []
        for mask in range(8):
            signs = [(-1)**((mask & j).bit_count()) for j in range(8)]
            expected.append(np.dot(x, signs)/(8 if mask == 0 else 4))
        got = fwt(x).real/np.array([8]+[4]*7)
        np.testing.assert_array_equal(got, expected)
        np.testing.assert_array_equal(got, [73.25, -23, 5, 1.5, -1.5, 10, 0, -.5])

    def test_macwilliams_prime_and_extension_fields(self):
        import galois
        X, Y = sp.symbols('X Y')
        for q in [2, 3, 4, 5, 8, 9]:
            field = galois.GF(q, compile='python-calculate')
            v = field([1, 1, q-1])
            code = [a*v for a in field.elements]
            dual = [field(t) for t in product(range(q), repeat=3)
                    if np.sum(field(t)*v) == 0]
            weight = lambda h: sum(int(a != 0) for a in h)
            polynomial = lambda rows: sum(X**(3-weight(h))*Y**weight(h) for h in rows)
            lhs = polynomial(dual)
            rhs = polynomial(code).subs({X: X+(q-1)*Y, Y: X-Y}, simultaneous=True)/q
            self.assertEqual(sp.expand(lhs-rhs), 0)

    def test_walsh_image_projection_against_basis_sum(self):
        rng = np.random.default_rng(4)
        x = rng.normal(size=(4, 8)) + 1j*rng.normal(size=(4, 8))
        rows = lambda n: np.array([[(-1)**((i&j).bit_count())
                                   for j in range(n)] for i in range(n)])
        Wy, Wx = rows(4), rows(8)
        sy = np.sum(Wy[:, 1:] != Wy[:, :-1], axis=1)
        sx = np.sum(Wx[:, 1:] != Wx[:, :-1], axis=1)
        order = sorted(product(range(4), range(8)), key=lambda t: (sy[t[0]]+sx[t[1]], sy[t[0]]))
        for M in [0, 1, 3, 15, 32]:
            expected = np.zeros_like(x)
            for i, j in order[:M]:
                basis = Wy[i, :, None]*Wx[j, None, :]
                expected += np.sum(x*basis)/x.size*basis
            np.testing.assert_allclose(walsh_approx2d(x, M), expected, atol=2e-15)

    def test_theta_identity_and_truncation_bound(self):
        import mpmath as mp
        with mp.workdps(60):
            for t in map(mp.mpf, ['0.1', '0.3', '1', '2', '10']):
                theta = lambda t: 1+2*sum(mp.exp(-mp.pi*n*n*t) for n in range(1, 101))
                self.assertLess(abs(theta(t)-theta(1/t)/mp.sqrt(t)), mp.mpf('1e-55'))
                if t < 1:
                    self.assertLess(theta(1/t)-1, 2*mp.exp(-mp.pi/t)/(1-mp.exp(-3*mp.pi/t)))


if __name__ == '__main__':
    unittest.main()
