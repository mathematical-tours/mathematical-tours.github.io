"""Independent checks of the corrected finite-group formulas.

These checks complement the proofs: they exercise sign and normalization
conventions, including cases for which the previous formulas fail.
"""
import unittest
from itertools import product
import numpy as np
from numpy.testing import assert_allclose


class FiniteGroupIdentities(unittest.TestCase):
    def test_character_extension_root(self):
        # Additive Z/NZ, H = d Z/NZ. A character of H sends d to
        # exp(2*pi*i*j/(N/d)); extending it across x=1 requires a d-th
        # root of this value, not an arbitrary d-th root of unity.
        for n in range(2,25):
            for d in range(2,n+1):
                if n%d:
                    continue
                for j in range(n//d):
                    zeta=np.exp(2j*np.pi*j/n)
                    chi_hd=np.exp(2j*np.pi*j/(n//d))
                    self.assertAlmostEqual(abs(zeta**d-chi_hd),0,places=12)
                    extension=np.array([np.exp(2j*np.pi*j*(a//d)/(n//d))*zeta**(a%d)
                                        for a in range(n)])
                    for a in range(n):
                        assert_allclose(extension[(a+np.arange(n))%n],
                                        extension[a]*extension,atol=3e-13)
        self.assertFalse(np.isclose(1**2,-1))  # Z/4Z counterexample.

    def test_fourier_convention_and_convolution(self):
        rng=np.random.default_rng(102)
        for n in [3,4,7,16,17]:
            x=np.arange(n)
            characters=np.exp(2j*np.pi*np.outer(x,x)/n)
            f=rng.normal(size=n)+1j*rng.normal(size=n)
            g=rng.normal(size=n)+1j*rng.normal(size=n)
            F=characters@f
            G=characters@g
            assert_allclose(characters.conj().T@F/n,f,atol=1e-12)
            assert_allclose(F,n*np.fft.ifft(f),atol=1e-12)
            assert_allclose(np.sum(f*g.conj()),np.sum(F*G.conj())/n,atol=1e-12)
            convolution=np.array([sum(f[h]*g[(k-h)%n] for h in range(n))
                                  for k in range(n)])
            assert_allclose(characters@convolution,F*G,atol=2e-11)
            assert_allclose(np.eye(n)@np.eye(n).T/n,np.eye(n)/n)

    def test_heisenberg_intertwining(self):
        n=7
        rng=np.random.default_rng(8)
        f=rng.normal(size=n)+1j*rng.normal(size=n)
        k=np.arange(n)
        F=np.exp(2j*np.pi*np.outer(k,k)/n)
        for x,j in product(range(n),repeat=2):
            lam=np.exp(.4j)
            chi=np.exp(2j*np.pi*j*k/n)
            action=lam*chi*f[(k+x)%n]
            expected=lam*np.exp(-2j*np.pi*j*x/n)*np.exp(-2j*np.pi*k*x/n)*(F@f)[(k+j)%n]
            assert_allclose(F@action,expected,atol=2e-13)

    def test_rotation_trace(self):
        for a in np.linspace(0,np.pi,31):
            c,s=np.cos(a),np.sin(a)
            r=np.array([[1,0,0],[0,c,-s],[0,s,c]])
            t=np.array([[c,0,s],[0,1,0],[-s,0,c]])
            assert_allclose((np.trace(r@t.T)-1)/2,((1+c)**2-2)/2,atol=5e-16)

    def test_indicator_bounds_and_uncertainty_extremizers(self):
        for n in [6,9,12]:
            k=np.arange(n)
            F=np.exp(2j*np.pi*np.outer(k,k)/n)
            for size in range(1,n//2+1):
                f=(k<size).astype(float)
                phi=np.max(abs(F@f)[1:])
                self.assertGreaterEqual(phi+1e-12,np.sqrt(size/2))
                self.assertLessEqual(phi,size+1e-12)
            for d in range(1,n+1):
                if n%d:
                    continue
                f=(k%d==0).astype(float)
                self.assertEqual(np.count_nonzero(f)*np.count_nonzero(abs(F@f)>1e-10),n)
                # Include both translation and modulation in an equality case.
                u=np.exp(2j*np.pi*k/n)*f[(k-1)%n]
                self.assertEqual(np.count_nonzero(abs(u)>1e-10)*
                                 np.count_nonzero(abs(F@u)>1e-10),n)


if __name__=='__main__':
    unittest.main()
