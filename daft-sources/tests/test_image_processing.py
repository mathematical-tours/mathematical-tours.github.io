import sys, tempfile, unittest
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.signal import convolve2d
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'listings'))
from image_processing import *


class ImageMethods(unittest.TestCase):
    def test_gaussian_and_reflection(self):
        kernel=gaussian_kernel(1.2)
        self.assertAlmostEqual(kernel.sum(),1)
        np.testing.assert_allclose(kernel,kernel[::-1,::-1])
        np.testing.assert_allclose(smooth_image(np.ones((12,16)),1.2),1,atol=1e-14)
        image=np.random.default_rng(1).normal(size=(17,19));r=kernel.shape[0]//2
        expected=convolve2d(np.pad(image,r,mode='reflect'),kernel,mode='valid')
        np.testing.assert_allclose(smooth_image(image,1.2),expected,atol=1e-14)
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'sample.png';Image.fromarray(np.full((11,13),128,dtype=np.uint8)).save(path)
            source,filtered=apply_gaussian(path)
            np.testing.assert_allclose(source,128/255);np.testing.assert_allclose(source,filtered)

    def test_sliding_sums_and_correlation(self):
        rng=np.random.default_rng(4);f=rng.normal(size=(15,18));g=f[7:11,9:13].copy()
        for power in [1,2]:
            expected=np.array([[np.sum(f[i:i+4,j:j+4]**power) for j in range(15)] for i in range(12)])
            np.testing.assert_allclose(sliding_sum(f,4,power),expected,atol=1e-13)
        rho=normalized_correlation(f,g)
        expected=np.empty_like(rho)
        for i,j in np.ndindex(rho.shape):
            a=f[i:i+4,j:j+4].ravel();b=g.ravel()
            expected[i,j]=np.corrcoef(a,b)[0,1]
        np.testing.assert_allclose(rho,expected,atol=1e-13)
        self.assertEqual(np.unravel_index(np.argmax(rho),rho.shape),(7,9))
        np.testing.assert_allclose(normalized_correlation(3*f+100,5*g-20),rho,atol=1e-12)
        with self.assertRaises(ValueError):normalized_correlation(f,np.ones((4,4)))
        self.assertTrue(np.all(np.isnan(normalized_correlation(np.ones((8,8)),g))))

    def test_fourier_resampling(self):
        rng=np.random.default_rng(3);x=rng.normal(size=(15,17))
        np.testing.assert_allclose(fft_translate(x,(2,-3)),np.roll(x,(2,-3),axis=(0,1)),atol=1e-13)
        np.testing.assert_allclose(fft_translate(fft_translate(x,(.3,1.2)),(-.3,-1.2)),x,atol=1e-13)
        for axis in [0,1]:
            y=fft_shear(x,.27,axis)
            self.assertAlmostEqual(np.linalg.norm(y),np.linalg.norm(x))
            np.testing.assert_allclose(fft_shear(y,-.27,axis),x,atol=1e-13)
        np.testing.assert_allclose(fft_rotate(fft_rotate(x,.3),-.3),x,atol=1e-13)
        error=np.linalg.norm(fft_rotate(fft_rotate(x,.2),.3)-fft_rotate(x,.5))
        self.assertGreater(error,1)
        self.assertEqual(rotation_examples(x).shape[0],4)

    def test_frequency_response_samples(self):
        n=65
        for epsilon in [0,.25,1]:
            kernel=smooth_lowpass_kernel(n,epsilon)
            response=np.fft.fft(kernel)
            self.assertLess(np.max(abs(response.imag)),1e-13)
            self.assertAlmostEqual(response[0].real,1)
            self.assertTrue(np.all(response.real>=-1e-13))
            self.assertTrue(np.all(response.real<=1+1e-13))
            omega,curve=filter_response(kernel,n*8)
            k=np.rint(np.fft.fftfreq(n)*n).astype(int)
            direct=np.exp(-1j*omega[:,None]*k[None,:])@kernel
            np.testing.assert_allclose(curve,direct,atol=1e-13)
        np.testing.assert_allclose(lowpass_kernel(n),smooth_lowpass_kernel(n,0))

if __name__=='__main__':unittest.main()
