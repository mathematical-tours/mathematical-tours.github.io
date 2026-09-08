"""Scientific invariants used by multiple interactive figures."""
import unittest
from numerics import *
from experiments_images import hadamard,inpaint_case
from experiments_learning import cert_case,clustering_case

class ScientificChecks(unittest.TestCase):
    def test_bird_overview_preserves_peaks_and_detail_samples(self):
        from experiments_signal import acquisition
        fs,y=C.bird()
        for zoom in [1,16]:
            result=acquisition('examples',{'source':'Birdsong','zoom':zoom})
            overview,detail=result['panels'];envelope=overview['series'][0];samples=detail['series'][0]
            self.assertAlmostEqual(min(envelope['y']),float(y.min()),places=7)
            self.assertAlmostEqual(max(envelope['y']),float(y.max()),places=7)
            start=int(len(y)*.4);count=len(samples['y'])
            np.testing.assert_allclose(samples['y'],y[start:start+count],atol=1e-8)
            np.testing.assert_allclose(samples['x'],np.arange(start,start+count)/fs,atol=1e-6)
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
            pts,weights,history,E=clustering_case(6,continuous)
            self.assertEqual(len(E),len(history))
            self.assertLessEqual(max(np.diff(E)),1e-10)
            for k in [0,1,2,3,40]:
                distances=np.sum((pts[:,None]-history[k][None])**2,axis=2).min(axis=1)
                self.assertAlmostEqual(E[k],np.average(distances,weights=weights),places=12)
    def test_empirical_spectrum_retains_all_eigenvalues(self):
        from experiments_learning import random_case
        for beta in [.1,.6,.8]:
            for P in [40,80]:
                panel=random_case('marcenko-pastur',{'beta':beta,'P':P})['panels'][0]
                hist=panel['series'][0]
                self.assertAlmostEqual(sum(hist['y'])*hist['barWidth'],1,places=7)
    def test_logistic_upper_bound_touches_binary_loss(self):
        from experiments_learning import learning
        panel=learning('losses',{'margin':0})['panels'][0]
        values={s['label']:np.array(s['y']) for s in panel['series']}
        self.assertTrue(np.all(values['Logistic / log 2']>=values['Binary']))
        self.assertEqual(values['Logistic / log 2'][300],1)
        self.assertLess(values['Logistic'][300],1)
    def test_curvelet_grid_is_tangent_to_frequency_direction(self):
        from experiments_signal import directional
        for angle in [0,30,90]:
            result=directional('curvelet-discretization',{'scale':16,'angle':angle})
            s=result['panels'][0]['series'][0]
            x=np.array(s['x']).reshape(13,7);y=np.array(s['y']).reshape(13,7)
            tangent=np.array([x[0,1]-x[0,0],y[0,1]-y[0,0]])
            from curvelet_frame import windows
            _,info=windows(128)
            actual=min(info[1:],key=lambda v:abs(np.log2(v['radius']/16))*10+abs(np.angle(np.exp(1j*(v['angle']-np.deg2rad(angle))))))
            normal=np.array([np.cos(actual['angle']),np.sin(actual['angle'])])
            self.assertAlmostEqual(tangent@normal,0,places=7)
            self.assertAlmostEqual(np.linalg.norm(tangent),.25,places=7)
    def test_small_numbers_survive_export(self):
        self.assertEqual(arr(np.array([1e-15,1e-9,0])),[1e-15,1e-9,0])
    def test_heat_preserves_mean_and_reduces_dirichlet_energy(self):
        f=input_image('Flower');g=C.heat(f,2)
        self.assertAlmostEqual(f.mean(),g.mean(),places=12)
        self.assertLess(np.sum(C.gradient(g)**2),np.sum(C.gradient(f)**2))

if __name__=='__main__':unittest.main()
