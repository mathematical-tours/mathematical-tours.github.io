import sys, itertools, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'listings'))
from splines import *


class ChapterFourMathematics(unittest.TestCase):
    def test_polygon_plancherel_identities(self):
        rng=np.random.default_rng(71)
        for n in [3,4,9,16]:
            z=rng.normal(size=n)+1j*rng.normal(size=n);F=np.fft.fft(z);k=np.arange(1,n)
            area=np.sum((z.conj()*np.roll(z,-1)).imag)/2
            T=np.sum(abs(z-z.mean())**2);S=np.sum(abs(z-np.roll(z,-1))**2)
            self.assertAlmostEqual(T,np.sum(abs(F[1:])**2)/n)
            self.assertAlmostEqual(area,np.sum(np.sin(2*np.pi*k/n)*abs(F[1:])**2)/(2*n))
            self.assertAlmostEqual(S,4*np.sum(np.sin(np.pi*k/n)**2*abs(F[1:])**2)/n)
            self.assertLessEqual(abs(area),T/2+1e-12)
            self.assertLessEqual(4*np.sin(np.pi/n)**2*T,S+1e-12)

    def test_spline_symbols_and_cardinal_values(self):
        omega=np.linspace(-np.pi,np.pi,1001)
        for degree in range(1,14):
            symbol=spline_symbol(degree,omega)
            self.assertTrue(np.all(symbol>0))
            self.assertTrue(np.all(symbol+1e-14>=np.sinc(omega/(2*np.pi))**(degree+1)))
        np.testing.assert_allclose(centered_bspline(3,np.arange(-2,3)),[0,1/6,2/3,1/6,0])
        np.testing.assert_allclose(spline_symbol(3,omega),(2+np.cos(omega))/3,atol=1e-14)
        k=np.arange(-12,13);np.testing.assert_allclose(cardinal_cubic(k),(k==0).astype(float),atol=1e-15)

    def test_cardinal_spline_limit(self):
        eta=np.linspace(-np.pi,np.pi,1501)[1:-1]
        frequencies=eta[None,:]+2*np.pi*np.arange(-8,9)[:,None]
        target=(abs(frequencies)<np.pi).astype(float)
        errors=[]
        for degree in [1,3,7,15]:
            H=cardinal_spline_spectrum(degree,frequencies)
            errors.append(np.mean((H-target)**2))
        self.assertTrue(np.all(np.diff(errors)<0))
        self.assertLess(errors[-1],errors[0]/4)

    def test_aliasing_identity_and_trapezoidal_accuracy(self):
        n=16;k=np.arange(n);f=np.exp(np.cos(2*np.pi*k/n))
        from scipy.special import iv
        coefficients=np.fft.fft(f)/n
        for j in [0,1,3,7]:
            alias=sum(iv(abs(j+ell*n),1) for ell in range(-8,9))
            self.assertAlmostEqual(coefficients[j].real,alias,places=14)
        self.assertLess(abs(coefficients[0]-iv(0,1)),1e-13)

    def test_candy_consensus_and_even_cycle(self):
        for n in [3,4,5]:
            for initial in itertools.product([0,2,4],repeat=n):
                x=np.array(initial)
                for _ in range(3*n):
                    if np.ptp(x)==0:break
                    x=2*((x+np.roll(x,1)+3)//4)
                self.assertEqual(np.ptp(x),0)
        x=np.array([0,2,0,2]);y=2*((np.roll(x,1)+np.roll(x,-1)+3)//4)
        self.assertFalse(np.array_equal(x,y))
        np.testing.assert_array_equal(2*((np.roll(y,1)+np.roll(y,-1)+3)//4),x)

if __name__=='__main__':unittest.main()
