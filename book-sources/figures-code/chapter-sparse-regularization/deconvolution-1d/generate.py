from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


t,f=signal_1d(512);H=blur_symbol(f.shape,2);h=fft.ifft(H).real;y=convolution(f,H)+.02*np.random.default_rng(20260906).standard_normal(f.shape);g=sparse_deblur(y,H,.005,steps=1200);fig,a=canvas(1,4,height=3.1,width=13)
for ax,z,title in zip(a.ravel(),[f,fft.fftshift(h),y,g],['Clean signal','Periodic Gaussian filter','Blurred noisy observation',f'Wavelet synthesis recovery: {snr(f,g):.1f} dB']):curve(ax,t,z,title,'Position','Amplitude')
finish(fig,__file__,'One-dimensional deconvolution solves .5||H Psi c-y||²+.005||c||1 in a periodic orthonormal db2 basis. The synthesis coefficients, not image values, are soft-thresholded by the optimizer. The same filter produces the observation and the forward/adjoint operators.',parameters={'N':512,'blur_sigma':2,'noise_sigma':.02,'seed':20260906,'lambda':.005,'iterations':1200})
