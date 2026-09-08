"""Scientific invariants used by multiple interactive figures."""
import unittest
from numerics import *
from experiments_images import hadamard,inpaint_case
from experiments_learning import cert_case,clustering_case

class ScientificChecks(unittest.TestCase):
    def test_gradient_divergence_adjoint(self):
        rng=np.random.default_rng(19);f=rng.normal(size=(17,19));p=rng.normal(size=(2,17,19))
        self.assertLess(abs(np.sum(C.gradient(f)*p)+np.sum(f*C.divergence(p))),1e-10)
    def test_wavelet_perfect_reconstruction_and_energy(self):
        f=input_image('Flower',128)
        for w in ['haar','db2','db4']:
            for level in [1,2,3,4]:
                a,s=C.wavelet_array(f,w,level);g=C.wavelet_inverse(a,s,w)
                self.assertLess(np.linalg.norm(g-f),1e-9);self.assertAlmostEqual(np.sum(a*a),np.sum(f*f),places=8)
    def test_equal_budget_orthonormal_approximation(self):
        for name in ['Flower','Cartoon']:
            for basis in ['Fourier','DCT','haar','db4']:
                a,inv,order=transform_image(name,basis);f=input_image(name)
                for M in [16,256,4096,16384]:
                    g,b=approximation_image(name,basis,M);linear,_=approximation_image(name,basis,M,False)
                    self.assertLessEqual(np.linalg.norm(f-g),np.linalg.norm(f-linear)+1e-10)
                    self.assertAlmostEqual(np.sum((f-g)**2),np.sum((a-b)**2),places=7)
    def test_hadamard_is_self_adjoint_unitary(self):
        rng=np.random.default_rng(21);x=rng.normal(size=256);y=rng.normal(size=256)
        self.assertLess(np.linalg.norm(hadamard(hadamard(x))-x),1e-12)
        self.assertAlmostEqual(hadamard(x)@y,x@hadamard(y),places=10)
    def test_inpainting_preserves_observed_samples(self):
        for method in ['Sobolev','TV']:
            for mask in ['Square','Stripes','Random']:
                f,m,y,g=inpaint_case(80,method,mask);self.assertLess(np.max(abs((g-f)[m])),1e-12)
    def test_huffman_bounds_and_prefix_property(self):
        for K in [2,3,4,8]:
            for a in [0,.5,2,4]:
                p=np.exp(-a*np.arange(K));p/=p.sum();codes=C.huffman_codes(dict(enumerate(p)));words=list(codes.values());H=-sum(p*np.log2(p));L=sum(p[k]*len(c) for k,c in codes.items())
                self.assertLessEqual(H,L+1e-10);self.assertLess(L,H+1+1e-10)
                self.assertTrue(all(not a.startswith(b) for a in words for b in words if a!=b))
    def test_tangent_aligned_meshes(self):
        from anisotropic_models import fitted_mesh,mesh_checks
        for k in [32,40,48,56,64]:
            p,t=fitted_mesh(k);c=mesh_checks(p,t)
            self.assertLess(c['jump_triangle_max_tangent_angle_degrees'],8)
            self.assertGreater(c['jump_triangle_min_aspect_ratio'],12)
            self.assertAlmostEqual(c['domain_area'],1,places=12)
    def test_fuchs_interpolates_support_signs(self):
        for s in [2,5,12]:
            A,x,y,I,J,eta=cert_case(s,24);self.assertLess(max(abs(eta[I]-np.sign(x[I]))),1e-9)
    def test_lloyd_energy_does_not_increase(self):
        for continuous in [False,True]:
            _,_,_,E=clustering_case(6,continuous);self.assertLessEqual(max(np.diff(E)),1e-10)
    def test_small_numbers_survive_export(self):
        self.assertEqual(arr(np.array([1e-15,1e-9,0])),[1e-15,1e-9,0])
    def test_heat_preserves_mean_and_reduces_dirichlet_energy(self):
        f=input_image('Flower');g=C.heat(f,2)
        self.assertAlmostEqual(f.mean(),g.mean(),places=12)
        self.assertLess(np.sum(C.gradient(g)**2),np.sum(C.gradient(f)**2))

if __name__=='__main__':unittest.main()
