import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'listings'))
from polynomials import *


class PolynomialArithmetic(unittest.TestCase):
    def test_karatsuba_exact_products(self):
        rng=np.random.default_rng(5)
        for m,n in [(1,1),(1,9),(3,5),(8,8),(19,26),(32,33)]:
            p=[int(v)*10**30 for v in rng.integers(-5,6,size=m)]
            q=[int(v) for v in rng.integers(-5,6,size=n)]
            expected=[sum(p[i]*q[k-i] for i in range(m) if 0<=k-i<n) for k in range(m+n-1)]
            self.assertEqual(karatsuba(p,q),expected)
        self.assertEqual(polynomial_add([1,2],[3,4,5]),[4,6,5])

    def test_base_conversion_and_fft_product(self):
        for base in [2,10,20,256]:
            for x,y in [(0,0),(0,123),(1,2),(10**90+17,10**70-1),(987654321,1234567)]:
                a,b=integer_to_digits(x,base),integer_to_digits(y,base)
                self.assertEqual(digits_to_integer(a,base),x)
                self.assertEqual(digits_to_integer(multiply_digits(a,b,base),base),x*y)
        self.assertTrue(integer_product_example())
        with self.assertRaises(ValueError):integer_to_digits(-1,10)
        with self.assertRaises(ValueError):digits_to_integer([10],10)

if __name__=='__main__':unittest.main()
