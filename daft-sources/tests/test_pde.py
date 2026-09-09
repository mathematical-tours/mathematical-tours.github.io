import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'listings'))
from pde import *


class PDESolvers(unittest.TestCase):
    def test_heat_stencils_and_stability(self):
        n=12;eye=np.eye(n);A=np.roll(eye,1,axis=0)+np.roll(eye,-1,axis=0)-2*eye
        x=np.random.default_rng(9).normal(size=n)
        np.testing.assert_allclose(heat_explicit(x,.2,7),np.linalg.matrix_power(eye+.2*A,7)@x)
        for theta in [0,.25,.5,1]:
            M=np.linalg.solve(eye-theta*.7*A,eye+(1-theta)*.7*A)
            np.testing.assert_allclose(heat_theta(x,.7,theta,4),np.linalg.matrix_power(M,4)@x,atol=1e-12)
        alternating=(-1.)**np.arange(n)
        self.assertGreater(np.linalg.norm(heat_theta(alternating,2,.25,8)),np.linalg.norm(alternating))
        for theta in [.5,1]:
            self.assertLessEqual(np.linalg.norm(heat_theta(alternating,50,theta,8)),np.linalg.norm(alternating))
        X=alternating[:,None]*alternating[None,:]
        np.testing.assert_allclose(heat_theta_2d(X,.25,0,3),-X,atol=1e-12)
        self.assertGreater(np.linalg.norm(heat_theta_2d(X,.3,0,3)),np.linalg.norm(X))

    def test_spectral_heat_modes_and_mean(self):
        n=65;t=np.arange(n)/n;x=2+np.cos(6*np.pi*t)+.2j*np.sin(4*np.pi*t)
        expected=2+np.exp(-.5*.01*(6*np.pi)**2)*np.cos(6*np.pi*t)+.2j*np.exp(-.5*.01*(4*np.pi)**2)*np.sin(4*np.pi*t)
        np.testing.assert_allclose(heat_periodic(x,.01),expected,atol=1e-13)
        grid,solutions=heat_example()
        np.testing.assert_allclose(solutions.mean(axis=1),heat_initial(grid).mean())
        self.assertTrue(np.all(np.diff(np.linalg.norm(solutions,axis=1))<=0))

    def test_sine_transform_and_odd_extensions(self):
        rng=np.random.default_rng(1)
        for n in [2,3,8]:
            index=np.arange(1,n);V=np.sin(np.pi*index[:,None]*index[None,:]/n)
            x=rng.normal(size=n-1)+1j*rng.normal(size=n-1)
            np.testing.assert_allclose(sine_transform(x),V@x,atol=1e-12)
            np.testing.assert_allclose(sine_transform(sine_transform(x))*2/n,x,atol=1e-12)
            X=rng.normal(size=(n-1,n-1));np.testing.assert_allclose(sine_transform_2d(X),V@X@V,atol=1e-12)
        x=rng.normal(size=(2,3,4));extended=odd_extend(x)
        for axis,n in enumerate(extended.shape):
            np.testing.assert_allclose(np.take(extended,(-np.arange(n))%n,axis=axis),-extended)

    def test_poisson_manufactured_solutions(self):
        for dimension in [2,3]:
            n=8;axes=np.meshgrid(*([np.linspace(0,1,n+1)]*dimension),indexing='ij')
            exact=sum(x*x for x in axes)+np.prod(axes,axis=0)
            rhs=np.full(exact.shape,2*dimension)
            calculated=poisson_dirichlet(rhs,exact,1/n)
            np.testing.assert_allclose(calculated,exact,atol=2e-13)
        errors=[]
        for n in [8,16,32]:
            _,solution,exact=poisson_example(n);errors.append(np.max(abs(solution-exact)))
        self.assertTrue(all(3.5<a/b<4.5 for a,b in zip(errors,errors[1:])))
        np.testing.assert_allclose(poisson_left_boundary(np.arange(5)),1)
        self.assertEqual(poisson_source(0,0),0)

    def test_dirichlet_eigenmode(self):
        n=11;h=1/n;i=np.arange(1,n)
        mode=np.sin(2*np.pi*i/n)[:,None]*np.sin(3*np.pi*i/n)[None,:]
        eigenvalue=-4/h**2*(np.sin(np.pi/n)**2+np.sin(3*np.pi/(2*n))**2)
        np.testing.assert_allclose(poisson_zero_boundary(eigenvalue*mode,h),mode,atol=1e-13)

if __name__=='__main__':unittest.main()
