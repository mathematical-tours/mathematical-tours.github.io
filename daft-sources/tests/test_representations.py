import sys,itertools,math,unittest
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'listings'))
from representations import finite_wavelet,inverse_finite_wavelet

class RepresentationChecks(unittest.TestCase):
    def test_affine_action_and_character(self):
        for p in [2,3,5,7]:
            def action(b,a):return np.eye(p)[(pow(a,-1,p)*(np.arange(p)-b))%p]
            group=list(itertools.product(range(p),range(1,p)))
            for b,a in group:
                for d,c in group:
                    np.testing.assert_array_equal(action(b,a)@action(d,c),action((b+a*d)%p,a*c%p))
            character=np.array([np.trace(action(b,a))-1 for b,a in group])
            self.assertEqual(sum(character**2),len(group))

    def test_wavelet_complex_direct_and_inverse(self):
        rng=np.random.default_rng(714)
        for p in [2,3,7,17,53]:
            f=rng.normal(size=p)+1j*rng.normal(size=p)
            psi=rng.normal(size=p)+1j*rng.normal(size=p);psi-=psi.mean()
            W=finite_wavelet(f,psi)
            direct=np.array([[np.sum(f*np.conj(psi[(pow(a,-1,p)*(np.arange(p)-b))%p]))
                             for b in range(p)] for a in range(1,p)])
            np.testing.assert_allclose(W,direct,atol=2e-12)
            np.testing.assert_allclose(inverse_finite_wavelet(W,psi),f-f.mean(),atol=3e-13)
            c=p*np.sum(abs(psi)**2)
            self.assertAlmostEqual(np.sum(abs(W)**2)/np.sum(abs(f-f.mean())**2),c,places=8)
            P=np.fft.fft(psi);P[0]=np.sqrt(np.sum(abs(P[1:])**2)/(p-1))*np.exp(.7j)
            tight=np.fft.ifft(P);W=finite_wavelet(f,tight)
            np.testing.assert_allclose(inverse_finite_wavelet(W,tight,'tight'),f,atol=3e-13)
            self.assertAlmostEqual(np.mean(abs(W)**2)/np.mean(abs(f)**2),c/(p-1),places=8)
            delta=np.eye(p)[0]
            np.testing.assert_allclose(inverse_finite_wavelet(finite_wavelet(f,delta),delta,'tight'),f,atol=3e-13)
        with self.assertRaises(ValueError):finite_wavelet(np.ones(4),np.ones(4))
        with self.assertRaises(ValueError):inverse_finite_wavelet(np.zeros((2,3)),np.ones(3))

    def test_hadamard_distance_and_formal_dual(self):
        words=['00000000000','11011100010','01101110001','10110111000','01011011100','00101101110',
               '00010110111','10001011011','11000101101','11100010110','01110001011','10111000101']
        for a,b in itertools.combinations(words,2):self.assertEqual(sum(x!=y for x,y in zip(a,b)),6)
        t=sp.Symbol('t');P=sp.Poly(sp.expand(((1+t)**11+11*(1+t)**5*(1-t)**6)/12),t)
        self.assertEqual([P.nth(i) for i in [0,3,4,5,6,7,8,11]],
                         [1,sp.Rational(55,3),sp.Rational(110,3),sp.Rational(88,3),sp.Rational(88,3),sp.Rational(110,3),sp.Rational(55,3),1])
        self.assertEqual(P.eval(1),sp.Rational(512,3))
        paley=['0000000','1001011','1100101','1110010','0111001','1011100','0101110','0010111']
        values={int(w,2) for w in paley}
        self.assertTrue(all(a^b in values for a,b in itertools.product(values,repeat=2)))

    def test_linear_programming_bound(self):
        def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0
        for n,d,expected in [(7,3,16),(8,4,16),(11,6,12)]:
            matrix=np.array([[sum((-1)**j*choose(i,j)*choose(n-i,k-j) for j in range(k+1))
                              for i in range(n+1)] for k in range(n+1)])
            np.testing.assert_array_equal(matrix.sum(axis=0),[2**n]+[0]*n)
            bounds=[(1,1)]+[(0,0)]*(d-1)+[(0,None)]*(n-d+1)
            result=linprog(-np.ones(n+1),A_ub=-matrix,b_ub=np.zeros(n+1),bounds=bounds,method='highs')
            self.assertTrue(result.success);self.assertAlmostEqual(-result.fun,expected)

    def test_necklace_count(self):
        words=list(itertools.product(range(2),repeat=6))
        orbits=set()
        for word in words:
            transforms=[tuple(word[(s*j+b)%6] for j in range(6)) for s in [-1,1] for b in range(6)]
            orbits.add(min(transforms))
        self.assertEqual(len(orbits),13)
        counts=[sum(all(word[j]==word[(s*j+b)%6] for j in range(6)) for word in words)
                for s in [-1,1] for b in range(6)]
        self.assertEqual(sum(counts),156)

if __name__=='__main__':unittest.main()
