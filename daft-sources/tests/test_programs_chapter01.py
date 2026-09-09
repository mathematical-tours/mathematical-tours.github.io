"""Execute every Python listing accompanying Chapter 1."""
from pathlib import Path
import sys
import unittest
import numpy as np
from numpy.testing import assert_allclose
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from listings.finite_groups import circulant_determinant, rotation_product_trace, random_walk


class PythonPrograms(unittest.TestCase):
    def test_circulant_determinant(self):
        rng=np.random.default_rng(41)
        for n in [1,2,5,10,17]:
            for complex_data in [False,True]:
                a=rng.normal(size=n)
                if complex_data:
                    a=a+1j*rng.normal(size=n)
                i=np.arange(n)
                matrix=a[(i[None,:]-i[:,None])%n]
                assert_allclose(circulant_determinant(a),np.linalg.det(matrix),rtol=3e-12,atol=1e-10)
        assert_allclose(circulant_determinant(np.ones(4)),0,atol=1e-14)

    def test_symbolic_rotation_trace(self):
        result=rotation_product_trace()
        a=sp.symbols('alpha',real=True)
        self.assertEqual(sp.simplify(result-(2*sp.cos(a)+sp.cos(a)**2)),0)
        self.assertEqual(result.subs(a,0),3)
        self.assertEqual(result.subs(a,sp.pi),-1)

    def test_random_walk(self):
        rng=np.random.default_rng(44)
        for n in [3,8,15]:
            p0=rng.random(n); p0/=p0.sum()
            v=rng.random(n); v/=v.sum()
            P=np.array([[v[(j-i)%n] for i in range(n)] for j in range(n)])
            for k in [0,1,2,7,50]:
                actual=random_walk(p0,v,k)
                assert_allclose(actual,np.linalg.matrix_power(P,k)@p0,atol=2e-14)
                self.assertAlmostEqual(actual.sum(),1,places=12)
                self.assertGreaterEqual(actual.min(),-1e-14)

if __name__=='__main__':
    unittest.main()
