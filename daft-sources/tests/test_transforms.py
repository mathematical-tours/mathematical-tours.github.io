"""Direct-matrix and basis-integral references for elementary transforms."""
from pathlib import Path
import sys
import unittest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from listings.transforms import fwt, fwt2d, haar, haar_inverse, walsh_sequency, tensor_transform
from listings.finite_groups import paley_hadamard, paley_example


def walsh_matrix(n):
    return np.array([[(-1)**((i&j).bit_count()) for j in range(n)] for i in range(n)])


def haar_matrix(n):
    # Explicit normalized wavelets, evaluated on cell midpoints, in the
    # packed fine-to-coarse order. This independently fixes the detail sign.
    x=(np.arange(n)+.5)/n
    rows=[]
    for scale in range(n.bit_length()-2,-1,-1):
        for k in range(2**scale):
            t=2**scale*x-k
            row=np.where((0<=t)&(t<.5),1,0)-np.where((.5<=t)&(t<1),1,0)
            rows.append(row*np.sqrt(2**scale/n))
    rows.append(np.ones(n)/np.sqrt(n))
    return np.array(rows)


class FastTransforms(unittest.TestCase):
    def test_fwt_and_sequency(self):
        rng=np.random.default_rng(1)
        for n in [1,2,4,8,16,32,64]:
            W=walsh_matrix(n)
            x=rng.normal(size=n)+1j*rng.normal(size=n)
            assert_allclose(fwt(x),W@x,atol=2e-13)
            assert_allclose(fwt(fwt(x))/n,x,atol=2e-13)
            counts=np.sum(W[:,:-1]!=W[:,1:],axis=1)
            assert_array_equal(walsh_sequency(n),counts)
            assert_array_equal(np.sort(counts),np.arange(n))

    def test_rectangular_fwt2d(self):
        rng=np.random.default_rng(2)
        for n,m in [(1,1),(2,8),(16,4),(8,8)]:
            X=rng.normal(size=(n,m))+1j*rng.normal(size=(n,m))
            assert_allclose(fwt2d(X),walsh_matrix(n)@X@walsh_matrix(m).T,atol=2e-13)
            assert_allclose(fwt2d(fwt2d(X))/(n*m),X,atol=2e-13)

    def test_haar_basis_sign_normalization_and_inverse(self):
        rng=np.random.default_rng(3)
        for n in [1,2,4,8,16,32,64]:
            H=haar_matrix(n)
            assert_allclose(H@H.T,np.eye(n),atol=1e-14)
            x=rng.normal(size=n)+1j*rng.normal(size=n)
            y=haar(x)
            assert_allclose(y,H@x,atol=2e-13)
            assert_allclose(haar_inverse(y),x,atol=2e-13)
            assert_allclose(y[-1],np.sqrt(n)*x.mean(),atol=2e-13)
            assert_allclose(np.linalg.norm(x),np.linalg.norm(y),atol=2e-13)

    def test_tensor_product(self):
        rng=np.random.default_rng(4)
        for s in [2,3,4]:
            A=rng.normal(size=(s,s))+1j*rng.normal(size=(s,s))
            matrix=np.ones((1,1),complex)
            for depth in range(4):
                if depth:
                    matrix=np.kron(A,matrix)
                x=rng.normal(size=s**depth)+1j*rng.normal(size=s**depth)
                assert_allclose(tensor_transform(x,A),matrix@x,atol=2e-12)

    def test_paley(self):
        for p in [3,7,11,19,31]:
            H=paley_hadamard(p)
            assert_array_equal(H@H.T,(p+1)*np.eye(p+1,dtype=int))
            self.assertEqual(set(H.ravel()),{-1,1})
            Q=H[1:,1:]+np.eye(p,dtype=int)
            assert_array_equal(Q,-Q.T)
            assert_array_equal(Q@Q.T,p*np.eye(p,dtype=int)-np.ones((p,p),int))
        assert_array_equal(paley_example(),32*np.eye(32,dtype=int))

    def test_invalid_lengths(self):
        for method in [fwt,haar,haar_inverse]:
            for x in [[],[1,2,3],np.ones((2,2))]:
                with self.assertRaises(ValueError): method(x)
        with self.assertRaises(ValueError): tensor_transform(np.ones(6),np.eye(2))
        with self.assertRaises(ValueError): paley_hadamard(15)

if __name__=='__main__':
    unittest.main()
