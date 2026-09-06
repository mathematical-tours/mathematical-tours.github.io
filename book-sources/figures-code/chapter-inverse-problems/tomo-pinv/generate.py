from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f=phantom(128);F=fft.fft2(f,norm='ortho');fig,a=canvas(1,3,height=3.7);image(a[0,0],f,'Original image');records={}
for ax,K in zip(a[0,1:],[13,32]):mask=radial_mask(128,K);g=fft.ifft2(F*mask,norm='ortho').real;image(ax,g,f'{K} radial lines\n{snr(f,g):.1f} dB');records[str(K)]={'snr_db':snr(f,g),'measurement_fraction':float(mask.mean()),'observed_fourier_error':float(np.max(abs((fft.fft2(g,norm='ortho')-F)*mask)))}
finish(fig,__file__,'Implements the chapter’s partial-Fourier pseudoinverse model by setting unobserved Fourier entries to zero. The mask includes conjugate frequencies, so reconstruction is real. This is the exact pseudoinverse of that Cartesian restriction operator; it is not mislabeled filtered backprojection or the exact pseudoinverse of an interpolation-based discrete Radon matrix.',parameters={'size':128,'radial_line_counts':[13,32]},checks=records)
