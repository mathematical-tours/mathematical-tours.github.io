import sys, itertools, math, unittest
from pathlib import Path
import numpy as np
import sympy as sp
import galois
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'listings'))
from finite_fields import *
from error_correcting_codes import *


class ExactTransformsAndCodes(unittest.TestCase):
    def test_cyclotomic_search_with_independent_coefficients(self):
        for n in range(1,251):
            self.assertEqual(cyclotomic_coefficients(n),
                             [int(c) for c in reversed(cyclotomic_example(n).all_coeffs())])
        expected = [1,105,385,1365,1785,2805,3135,6545,6545,
                    10465,10465,10465,10465,10465,11305,11305,11305,11305,
                    11305,11305,11305,15015,11305,17255,17255,20615,20615]
        result = first_cyclotomic_coefficients(27,20615)
        self.assertEqual([result[a] for a in range(1,28)],expected)

    def test_cyclotomic_examples_and_factor_degrees(self):
        x=sp.Symbol('X')
        self.assertEqual(cyclotomic_example(6).as_expr(),x*x-x+1)
        self.assertIn(-2,cyclotomic_example(105).all_coeffs())
        self.assertEqual(first_cyclotomic_coefficients(2,105),{1:1,2:105})
        for p,n in [(2,15),(2,31),(3,8),(5,16),(7,9),(17,16)]:
            F,root=field_with_root(p,n)
            degree=int(sp.n_order(p,n))
            self.assertEqual(F.degree,degree)
            self.assertEqual(root**n,1)
            for ell in sp.factorint(n):self.assertNotEqual(root**(n//ell),1)
        with self.assertRaises(ValueError):field_with_root(2,16)

    def test_exact_finite_transforms(self):
        finite_transform_example()
        for p,n in [(2,15),(2,31),(3,8),(5,16),(7,9),(17,16)]:
            F,root=field_with_root(p,n);a=F(np.arange(n)%F.order);b=F((3*np.arange(n)+1)%F.order)
            transformed=finite_dft(a,root)
            np.testing.assert_array_equal(finite_dft(transformed,root,True),a)
            np.testing.assert_array_equal(finite_fft(a,root),transformed)
            convolution=F.Zeros(n)
            for k in range(n):
                for j in range(n):convolution[k]+=a[j]*b[(k-j)%n]
            np.testing.assert_array_equal(finite_dft(convolution,root),transformed*finite_dft(b,root))

    def test_principal_roots_and_modular_fft(self):
        for n,modulus,root in [(2,3,2),(4,5,2),(8,17,2),(16,257,2),(64,2**32+1,2)]:
            a=[j*j-2*j+3 for j in range(n)]
            direct=[sum(a[j]*pow(root,(-j*k)%n,modulus) for j in range(n))%modulus for k in range(n)]
            np.testing.assert_array_equal(modular_fft(a,modulus,root),direct)
            np.testing.assert_array_equal(modular_fft(direct,modulus,root,True),[v%modulus for v in a])
        with self.assertRaises(ValueError):modular_fft([1,2,3,4],15,2)
        # All residues and root orders for small composite moduli validate the criterion.
        for modulus in range(2,70):
            primes=list(sp.factorint(modulus))
            for n in range(1,9):
                candidates=[z for z in range(modulus) if pow(z,n,modulus)==1 and
                            all(math.gcd(1-pow(z,k,modulus),modulus)==1 for k in range(1,n))]
                self.assertEqual(bool(candidates),all((p-1)%n==0 for p in primes))

    def test_integer_products_without_rounding(self):
        for base in [2,10,10000]:
            for a,b in [(0,13),(1,1),(-17,29),(10**57+381,2**193-13),(-10**801+73,10**609-113)]:
                self.assertEqual(multiply_integers_modular(a,b,base),a*b)

    def test_reed_muller_and_boolean_nonlinearity(self):
        for n in [2,3,4]:
            for a in range(2**n):
                for b in [0,1]:
                    word=encode_rm(a,b,n)
                    direct=np.array([((a&x).bit_count()+b)%2 for x in range(2**n)])
                    np.testing.assert_array_equal(word,direct)
                    np.testing.assert_array_equal(decode_rm(word),word)
                    if n>=3:
                        for j in [0,2**n-1]:
                            noisy=word.copy();noisy[j]^=1
                            np.testing.assert_array_equal(decode_rm(noisy),word)
        # Including the zero-frequency coefficient is necessary for constant functions.
        self.assertEqual(4-np.max(abs(fwt(np.ones(8))))/2,0)
        for n in [2,4,6]:
            signs=np.array([(-1)**sum(((x>>(2*j))&1)*((x>>(2*j+1))&1) for j in range(n//2)) for x in range(2**n)])
            np.testing.assert_array_equal(abs(fwt(signs)),np.full(2**n,2**(n//2)))

    def test_hamming_decoding_and_distributions(self):
        for bits in [2,3,4]:
            H=hamming_check(bits);n=H.shape[1]
            if bits<=3:
                code=[np.array(word) for word in itertools.product([0,1],repeat=n) if np.all(H@word%2==0)]
                counts=np.bincount([sum(word) for word in code],minlength=n+1)
                np.testing.assert_array_equal(counts,hamming_weight_distribution(bits))
                for word in code:
                    for index in range(n):
                        noisy=word.copy();noisy[index]^=1
                        np.testing.assert_array_equal(decode_hamming(noisy),word)
            A=hamming_weight_distribution(bits)
            self.assertEqual(sum(A),2**(n-bits));self.assertEqual(A[3],n*(n-1)//6)
            self.assertEqual(A,A[::-1]);self.assertEqual(A[1:3],[0,0])
        np.testing.assert_array_equal(binary_digits(13,4),[1,0,1,1])

    def test_cyclic_hamming_and_extended_self_duality(self):
        G=np.array([[1,1,0,1,0,0,0],[0,1,1,0,1,0,0],[0,0,1,1,0,1,0],[0,0,0,1,1,0,1]]).T
        H=np.array([[1,0,0,1,0,1,1],[0,1,0,1,1,1,0],[0,0,1,0,1,1,1]])
        np.testing.assert_array_equal(H@G%2,0)
        code={(tuple(G@u%2)) for u in itertools.product([0,1],repeat=4)}
        self.assertEqual(len(code),16)
        for word in code:self.assertIn(tuple(np.roll(word,1)),code)
        extended=np.vstack([G,G.sum(axis=0)%2]);np.testing.assert_array_equal(extended.T@extended%2,0)
        weights=[sum(extended@u%2) for u in itertools.product([0,1],repeat=4)]
        np.testing.assert_array_equal(np.bincount(weights,minlength=9),[1,0,0,0,14,0,0,0,1])

    def test_bch_generator_table(self):
        F,alpha=field_with_root(2,31,[1,1,0,1,1,1])
        representatives=[0,1,3,5,11,14,15];product=galois.Poly.One(F)
        for index in representatives:
            orbit=cyclotomic_coset(2,31,index);P=galois.Poly.Roots(alpha**np.array(orbit));product*=P
            self.assertTrue(np.all(P.coeffs**2==P.coeffs))
            self.assertEqual(P.degree,len(orbit))
        np.testing.assert_array_equal(product.coeffs,F([1]+[0]*30+[1]))
        code=bch_code(2,31,5,[1,1,0,1,1,1])
        np.testing.assert_array_equal(code['generator'],[1,1,1,1,0,1,1,0,1,1,1])
        self.assertEqual(code['k'],21)

    def test_bch_all_error_counts_including_zero(self):
        rng=np.random.default_rng(632)
        for p,n,d in [(2,15,5),(2,31,7),(3,13,5),(5,8,5)]:
            code=bch_code(p,n,d);word=encode_bch(rng.integers(0,p,code['k']),code)
            for count in range(code['t']+1):
                for trial in range(8):
                    positions=rng.choice(n,count,replace=False)
                    received=word.copy();received[positions]+=type(word)(rng.integers(1,p,count))
                    np.testing.assert_array_equal(decode_bch(received,code),word)
                    np.testing.assert_array_equal(decode_bch(received,code,method='roots'),word)
        bch_example()

    def test_bch_exhaustive_binary_two_error_patterns(self):
        code=bch_code(2,15,5);word=encode_bch(np.arange(code['k'])%2,code)
        for count in [0,1,2]:
            for positions in itertools.combinations(range(15),count):
                received=word.copy();received[list(positions)]+=type(word)(1)
                np.testing.assert_array_equal(decode_bch(received,code),word)

if __name__=='__main__':unittest.main()
