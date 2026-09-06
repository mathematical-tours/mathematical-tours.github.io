from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(256);mask=np.outer(np.hanning(256),np.hanning(256));g=f*mask
F=fft.fftshift(fft.fft2(f,norm='ortho'));G=fft.fftshift(fft.fft2(g,norm='ortho'))
scale=max(np.abs(F).max(),np.abs(G).max());db=lambda v:20*np.log10(np.maximum(abs(v)/scale,1e-6))
fig,axs=canvas(1,4,height=3.2)
image(axs[0,0],f,'Original image');a=image(axs[0,1],db(F),'Fourier magnitude (dB)',vmin=-80,vmax=0,cmap='magma')
image(axs[0,2],g,'Hann-masked image');image(axs[0,3],db(G),'Masked spectrum (dB)',vmin=-80,vmax=0,cmap='magma')
fig.colorbar(a,ax=[axs[0,1],axs[0,3]],shrink=.65,ticks=[-80,-40,0])
finish(fig,__file__,'Both centered spectra are recomputed from the displayed flower images and share a decibel reference and color scale. A separable Hann mask forces the image toward zero on the boundary, reducing discontinuities in its periodic extension. Windowing also changes the image and broadens spectral features; it is not a lossless correction.',data_sources=['data/flower.png'],parameters={'mask':'separable Hann','spectrum':'unitary FFT, 20 log10 amplitude, shared reference','db_limits':[-80,0]})
