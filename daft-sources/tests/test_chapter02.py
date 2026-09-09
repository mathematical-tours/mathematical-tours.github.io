"""Checks for Gauss sums, finite-field duality, and quadratic reciprocity."""
import os
from pathlib import Path
os.environ.setdefault('NUMBA_CACHE_DIR',str(Path(__file__).resolve().parents[1]/'build/numba-cache'))
import unittest
import numpy as np
from numpy.testing import assert_allclose
import galois


class GaussSums(unittest.TestCase):
    def test_gauss_identities_over_extension_fields(self):
        for q in [3,4,5,7,8,9,16,25,27]:
            with self.subTest(q=q):
                F=galois.GF(q,compile='python-calculate')
                x=F.elements
                p=F.characteristic
                trace=(x[:,None]*x[None,:]).field_trace()
                additive=np.exp(2j*np.pi*np.asarray(trace,dtype=int)/p)
                assert_allclose(additive@additive.conj().T,q*np.eye(q),atol=2e-12)
                alpha=F.primitive_element
                powers=np.array([int(alpha**k) for k in range(q-1)])
                multiplicative=np.zeros((q-1,q),complex)
                multiplicative[:,powers]=np.exp(2j*np.pi*np.outer(np.arange(q-1),np.arange(q-1))/(q-1))
                G=multiplicative@additive.T
                assert_allclose(G[0,0],q-1,atol=2e-12)
                assert_allclose(G[0,1:],-1,atol=2e-12)
                assert_allclose(G[1:,0],0,atol=2e-12)
                assert_allclose(abs(G[1:,1:]),np.sqrt(q),atol=2e-12)
                # The zero extension has normalized norm squared (q-1)/q.
                assert_allclose(np.mean(abs(G)**2,axis=1),q-1,atol=3e-12)
                for j in range(1,q-1):
                    opposite=(-j)%(q-1)
                    assert_allclose(G[j,1:]*G[opposite,1:],q*multiplicative[j,int(-F(1))],atol=5e-12)
                if p!=2:
                    s=F.degree
                    expected=(-1)**(s-1)*1j**(s*(p-1)**2//4)*np.sqrt(q)
                    assert_allclose(G[(q-1)//2,1],expected,atol=2e-12)
                    if q==9:
                        assert_allclose(G[4,1],3,atol=1e-12)

    def test_determinant_phase(self):
        for p in [3,5,7,11,13,19,23]:
            k=np.arange(1,p)
            T=np.exp(2j*np.pi*np.outer(k,k)/p)
            phase,logabs=np.linalg.slogdet(T)
            expected_phase=(-1)**((p-1)//2)*1j**((p-1)*(p-2)//2)
            assert_allclose(phase,expected_phase,atol=2e-12)
            assert_allclose(logabs,(p-2)/2*np.log(p),atol=2e-12)

    def test_reciprocity(self):
        primes=[3,5,7,11,13,17,19,23,29,31,37,41,43,47]
        def legendre(a,p):
            a%=p
            return 0 if not a else (1 if a in {k*k%p for k in range(1,p)} else -1)
        for p in primes:
            for r in primes:
                if p!=r:
                    self.assertEqual(legendre(p,r)*legendre(r,p),(-1)**((p-1)*(r-1)//4))


if __name__=='__main__':
    unittest.main()
