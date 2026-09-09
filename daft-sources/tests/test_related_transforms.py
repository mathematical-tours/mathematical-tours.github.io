import sys, unittest
from pathlib import Path
import numpy as np
from scipy.linalg import toeplitz
from scipy.ndimage import spline_filter1d
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'listings'))
from hartley import *
from chirp import *
from splines import cubic_prefilter_parallel, cubic_prefilter_cascade


class RelatedTransforms(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(51)

    def test_hartley_kernel_and_inversion(self):
        for n in [1, 2, 3, 8, 17, 64]:
            angle=2*np.pi*np.outer(np.arange(n), np.arange(n))/n
            H=np.cos(angle)+np.sin(angle)
            for x in [self.rng.normal(size=n), self.rng.normal(size=n)+1j*self.rng.normal(size=n)]:
                np.testing.assert_allclose(hartley(x), H@x, atol=4e-13)
                np.testing.assert_allclose(hartley(hartley(x))/n, x, atol=1e-14)
                if n & (n-1) == 0:
                    np.testing.assert_allclose(fht(x), H@x, atol=4e-13)
        x=self.rng.normal(size=32);F=np.fft.fft(x)
        np.testing.assert_allclose(hartley(x), F.real-F.imag, atol=1e-14)

    def test_hartley_convolution_and_correlation(self):
        for n in [1, 3, 16, 31]:
            a=self.rng.normal(size=n);b=self.rng.normal(size=n)
            direct=np.array([sum(a[j]*b[(k-j)%n] for j in range(n)) for k in range(n)])
            np.testing.assert_allclose(hartley_convolution(a,b),direct,atol=1e-14)
            c=hartley(a);reverse=(-np.arange(n))%n
            np.testing.assert_allclose(hartley(hartley_convolution(a,a[reverse])),
                                       (c*c+c[reverse]**2)/2,atol=1e-12)
            np.testing.assert_allclose(hartley_convolution(a, np.zeros(n)),0,atol=1e-14)

    def test_unitary_hartley_powers_and_tensor(self):
        x=self.rng.normal(size=17)+1j*self.rng.normal(size=17)
        np.testing.assert_allclose(fractional_hartley(x,0),x,atol=1e-14)
        np.testing.assert_allclose(fractional_hartley(x,1),hartley(x)/np.sqrt(17),atol=1e-14)
        for a in [-1.2,.3,.5,2.5]:
            y=fractional_hartley(x,a)
            self.assertAlmostEqual(np.linalg.norm(y),np.linalg.norm(x))
            np.testing.assert_allclose(fractional_hartley(y,-a),x,atol=1e-14)
        x=self.rng.normal(size=(7,12))+1j*self.rng.normal(size=(7,12))
        np.testing.assert_allclose(hartley2d(hartley2d(x))/84,x,atol=1e-14)

    def test_generalized_hartley_inverse(self):
        for n in [3,4,17]:
            angle=2*np.pi*np.outer(np.arange(n),np.arange(n))/n
            for lam in [.1,.6,2.1]:
                H=np.cos(angle+lam); inverse=2*np.cos(angle+np.pi/2-lam)/(n*np.sin(2*lam))
                np.testing.assert_allclose(H@inverse,np.eye(n),atol=5e-14)

    def test_hartley_signed_interpolation(self):
        n,m=17,96;k=np.rint(np.fft.fftfreq(n)*n).astype(int)
        x=self.rng.normal(size=n);g=hartley(x);padded=np.zeros(m);padded[k%m]=g
        values=hartley(padded)/n
        direct=np.sum(g[:,None]*(np.cos(2*np.pi*k[:,None]*np.arange(m)/m)
                                +np.sin(2*np.pi*k[:,None]*np.arange(m)/m)),axis=0)/n
        np.testing.assert_allclose(values,direct,atol=1e-14)
        from spectral_methods import trigonometric_interpolate
        # Direct Fourier formula independently fixes the signed-frequency convention.
        expected=np.fft.fft(x) @ np.exp(2j*np.pi*k[:,None]*np.arange(m)/m)/n
        np.testing.assert_allclose(values,expected,atol=1e-14)

    def test_finite_field_hartley(self):
        import galois
        finite_hartley_example()
        F=galois.GF(17,compile='python-calculate');root=F.primitive_element;iota=root**4
        x=F(np.arange(16)%17);powers=root**np.outer(np.arange(16),np.arange(16))
        H=(powers+powers**-1)/F(2)+(powers-powers**-1)/(F(2)*iota)
        y=finite_fht(x,root,iota)
        np.testing.assert_array_equal(y,H@x)
        np.testing.assert_array_equal(finite_fht(y,root,iota),F(16)*x)

    def test_chirp_rectangular_and_complex_spirals(self):
        for n,m in [(1,1),(1,7),(7,1),(3,8),(13,13),(17,31),(64,64)]:
            x=self.rng.normal(size=n)+1j*self.rng.normal(size=n)
            for z in [1+0j, np.exp(.17j), np.exp(.001+.13j), np.exp(-.001+.31j)]:
                exact=z**(-np.outer(np.arange(m),np.arange(n)))@x
                np.testing.assert_allclose(czt(x,z,m),exact,rtol=2e-11,atol=2e-11)
            np.testing.assert_allclose(czt(x,np.exp(2j*np.pi/n)),np.fft.fft(x),atol=3e-12)

    def test_toeplitz_embedding(self):
        for m,n in [(1,1),(1,9),(6,1),(3,8),(8,3),(13,17)]:
            c=self.rng.normal(size=m)+1j*self.rng.normal(size=m)
            r=self.rng.normal(size=n)+1j*self.rng.normal(size=n);r[0]=c[0]
            x=self.rng.normal(size=n)+1j*self.rng.normal(size=n)
            np.testing.assert_allclose(toeplitz_product(c,r,x),toeplitz(c,r)@x,atol=1e-14)

    def test_rader_prime_lengths(self):
        from sympy import primitive_root
        for p in [2,3,5,7,17,31,127]:
            x=self.rng.normal(size=p)+1j*self.rng.normal(size=p)
            np.testing.assert_allclose(rader_fft(x,int(primitive_root(p))),np.fft.fft(x),atol=4e-14)
        self.assertEqual(inverse_modulo(3,17),6)
        with self.assertRaises(ValueError):rader_fft(np.ones(7),2)

    def test_scaled_fourier_and_zoom(self):
        x=self.rng.normal(size=31);n=x.size
        for alpha in [0,-1,.23,1,2]:
            expected=np.exp(-2j*np.pi*alpha*np.outer(np.arange(n),np.arange(n))/n)@x
            np.testing.assert_allclose(scaled_fourier(x,alpha),expected,atol=1e-12)
        np.testing.assert_allclose(scaled_fourier(x,0),x.sum(),atol=1e-13)
        x=self.rng.normal(size=(7,11))
        np.testing.assert_allclose(scaled_fourier2d(x,1),np.fft.fft2(x),atol=1e-13)
        nu=2.1+.003*np.arange(123)
        expected=np.exp(-2j*np.pi*np.outer(nu,np.arange(n))/n)@np.arange(n)
        np.testing.assert_allclose(zoom_spectrum(np.arange(n),2.1,.003,123),expected,atol=1e-10)

    def test_simpson_quadrature_and_physical_phase(self):
        a,h=-1.2,.04;t=a+h*np.arange(51);xi0,dxi=.31,.07;m=29
        samples=np.exp(t);xi=xi0+dxi*np.arange(m)
        weights=np.ones(51);weights[1:-1:2]=4;weights[2:-1:2]=2
        expected=h/3*(np.exp(-1j*np.outer(xi,t))@(weights*samples))
        actual=simpson_fourier(samples,a,h,xi0,dxi,m)
        np.testing.assert_allclose(actual,expected,atol=1e-13)
        exact=(np.exp((1-1j*xi)*t[-1])-np.exp((1-1j*xi)*a))/(1-1j*xi)
        self.assertLess(np.max(abs(actual-exact)),3e-6)
        # Simpson integrates cubics exactly at zero frequency, including both endpoints.
        val=simpson_fourier(t**3,a,h,0,0,1)[0]
        self.assertAlmostEqual(val.real,(t[-1]**4-a**4)/4)

    def test_mirror_cubic_prefilters(self):
        for n in [1,2,3,8,31,100]:
            for u in [np.ones(n),self.rng.normal(size=n), self.rng.normal(size=n)+1j*self.rng.normal(size=n)]:
                c=cubic_prefilter_parallel(u);d=cubic_prefilter_cascade(u)
                np.testing.assert_allclose(c,d,atol=3e-15)
                if n==1: reconstructed=c
                else:
                    padded=np.pad(c,(1,1),mode='reflect')
                    reconstructed=(padded[:-2]+4*padded[1:-1]+padded[2:])/6
                np.testing.assert_allclose(reconstructed,u,atol=2e-15)
                if np.isrealobj(u):
                    np.testing.assert_allclose(c,spline_filter1d(u,order=3,mode='mirror'),atol=4e-15)

    def test_spectral_shortening_identity(self):
        n=127;beta=2.23;b=int(np.floor(beta));L=n*b/beta;r=round(L);s=r-L
        x=np.exp(2j*np.pi*beta*np.arange(r)/n)
        result=scaled_fourier(x,r/L);ell=b-np.arange(r);nonzero=ell!=0
        expected=abs(np.sin(np.pi*s*ell[nonzero]/L)/np.sin(np.pi*ell[nonzero]/L))
        np.testing.assert_allclose(abs(result[nonzero]),expected,atol=2e-11)
        self.assertAlmostEqual(abs(result[b]),r,places=10)
        near=(abs(ell)>0)&(abs(ell)<=L/2)
        self.assertLessEqual(np.max(abs(result[near])),np.pi*abs(s)/2+1e-11)

if __name__=='__main__':unittest.main()
