"""Test the mathematical normalization, branch choices, and eigenbasis repair."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.fft import dct

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'listings'))
from spectral_methods import (trigonometric_interpolate,dct2,dct3,
    chebyshev_interpolate,chebyshev_example,fractional_derivative,
    fourier_components,fractional_dft,fourier_eigenbasis,orthogonalize_translates)


class SpectralMethods(unittest.TestCase):
    def test_trigonometric_interpolation_complex_and_real(self):
        rng=np.random.default_rng(88)
        for n in [1,3,7,13]:
            for eta in [1,2,5]:
                x=rng.normal(size=n)+1j*rng.normal(size=n)
                frequencies=np.arange(-(n//2),n//2+1)
                spectrum=np.fft.fft(x)
                expected=np.exp(2j*np.pi*np.outer(np.arange(n*eta)/(n*eta),frequencies))@spectrum[frequencies%n]/n
                np.testing.assert_allclose(trigonometric_interpolate(x,eta),expected,atol=4e-15)
                np.testing.assert_allclose(trigonometric_interpolate(x,eta)[::eta],x,atol=1e-15)
                self.assertTrue(np.isrealobj(trigonometric_interpolate(x.real,eta)))

    def test_cosine_transforms_against_defining_sums(self):
        rng=np.random.default_rng(31)
        for n in [1,2,5,16,31]:
            x=rng.normal(size=n)
            C=np.cos(np.pi*np.outer(np.arange(n),np.arange(n)+.5)/n)
            modified=x.copy();modified[0]/=2
            np.testing.assert_allclose(dct2(x),C@x,atol=5e-14)
            np.testing.assert_allclose(dct3(x),C.T@modified,atol=5e-14)
            np.testing.assert_allclose(2*dct3(dct2(x))/n,x,atol=1e-15)
            np.testing.assert_allclose(2*dct2(x),dct(x,type=2),atol=5e-15)
            np.testing.assert_allclose(2*dct3(x),dct(x,type=3),atol=5e-15)

    def test_chebyshev_polynomial_reproduction(self):
        nodes=np.cos((np.arange(13)+.5)*np.pi/13)
        f=lambda x:2*x**10-3*x**7+.6*x**2-1
        grid=np.linspace(-1,1,121)
        np.testing.assert_allclose(chebyshev_interpolate(f(nodes),grid),f(grid),atol=9e-15)
        nodes,values,grid,result=chebyshev_example()
        np.testing.assert_allclose(chebyshev_interpolate(values,nodes),values,atol=3e-14)
        self.assertTrue(np.all(np.isfinite(result)))

    def test_fractional_derivative_on_fourier_modes(self):
        n,period=101,7.3;t=np.arange(n)*period/n
        for k in [-9,-1,0,2,11]:
            mode=np.exp(2j*np.pi*k*t/period)
            for alpha in [0,.4,1,1.7,2]:
                multiplier=1 if alpha==0 else abs(2*np.pi*k/period)**alpha*np.exp(.5j*np.pi*alpha*np.sign(k))
                np.testing.assert_allclose(fractional_derivative(mode,alpha,period),multiplier*mode,atol=8e-12)
        real=np.sin(4*np.pi*t/period)+.3*np.cos(6*np.pi*t/period)
        derivative=4*np.pi/period*np.cos(4*np.pi*t/period)-.3*6*np.pi/period*np.sin(6*np.pi*t/period)
        np.testing.assert_allclose(fractional_derivative(real,1,period),derivative,atol=8e-14)

    def test_projectors_and_fractional_group_law(self):
        rng=np.random.default_rng(43)
        for n in [1,2,3,4,7,12,31]:
            x=rng.normal(size=n)+1j*rng.normal(size=n)
            components=fourier_components(x)
            np.testing.assert_allclose(components.sum(axis=0),x,atol=1e-15)
            for component,eigenvalue in zip(components,[1,-1,1j,-1j]):
                np.testing.assert_allclose(np.fft.fft(component)/np.sqrt(n),eigenvalue*component,atol=1e-15)
            gram=components.conj()@components.T
            np.testing.assert_allclose(gram,np.diag(np.diag(gram)),atol=6e-15)
            for turns in [(0,0,0,0),(1,-2,3,0)]:
                for alpha,beta in [(.3,.7),(-.7,1.3),(2.8,-.4)]:
                    a=fractional_dft(fractional_dft(x,alpha,turns),beta,turns)
                    b=fractional_dft(x,alpha+beta,turns)
                    np.testing.assert_allclose(a,b,atol=3e-14)
                    self.assertAlmostEqual(np.linalg.norm(a),np.linalg.norm(x),places=12)
            np.testing.assert_allclose(fractional_dft(x,1),np.fft.fft(x)/np.sqrt(n),atol=1e-15)
            np.testing.assert_allclose(fractional_dft(x,2),x[(-np.arange(n))%n],atol=2e-15)

    def test_parity_eigenbasis_and_repeated_eigenvalue(self):
        for n in range(3,41):
            V=fourier_eigenbasis(n)
            U=np.exp(-2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)/np.sqrt(n)
            F=V.T@U@V
            np.testing.assert_allclose(V.T@V,np.eye(n),atol=3e-15)
            np.testing.assert_allclose(F,np.diag(np.diag(F)),atol=4e-14)
            np.testing.assert_allclose(np.diag(F)**4,1,atol=4e-14)
        S=np.roll(np.eye(4),1,axis=0)+np.roll(np.eye(4),-1,axis=0)+np.diag([-2,-4,-6,-4])
        for v in [np.array([0,1,0,-1]),np.array([1,-1,-1,-1])]:
            np.testing.assert_array_equal(S@v,-4*v)

    def test_translate_orthogonalization_normalized_gram(self):
        rng=np.random.default_rng(93)
        for n in [3,8,17,32]:
            x=rng.normal(size=n)+1j*rng.normal(size=n)
            y=orthogonalize_translates(x)
            orbit=np.column_stack([np.roll(y,k) for k in range(n)])
            np.testing.assert_allclose(orbit.conj().T@orbit/n,np.eye(n),atol=9e-16)
            phi=rng.normal(size=n)+1j*rng.normal(size=n)
            coefficients=np.fft.ifft(np.fft.fft(phi)*np.fft.fft(y).conj())/n
            np.testing.assert_allclose(orbit@coefficients,phi,atol=4e-15)


if __name__=='__main__': unittest.main()
