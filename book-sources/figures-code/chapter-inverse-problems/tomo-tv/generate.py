from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f=phantom(128);mask=radial_mask(128,13);y=fft.fft2(f,norm='ortho')*mask;pinv=fft.ifft2(y,norm='ortho').real;g=fourier_tv(y,mask,.02,1000);fig,a=canvas(1,3,height=3.7)
for ax,z,title in zip(a.ravel(),[f,pinv,g],['Original',f'Pseudoinverse: {snr(f,pinv):.1f} dB',f'TV: {snr(f,g):.1f} dB']):image(ax,z,title)
E=lambda x:.5*np.sum(abs(fft.fft2(x,norm='ortho')*mask-y)**2)+.02*tv(x)
assert E(g)<E(pinv)
finish(fig,__file__,'Minimizes .5||mask F x-y||²+.02 TV(x) for the same 13-line Cartesian Fourier acquisition as the pseudoinverse. The primal-dual method uses an exact Fourier-domain proximal step for fidelity and a pointwise Euclidean dual-ball projection for isotropic TV. The objective is checked against the initial pseudoinverse.',parameters={'size':128,'lines':13,'lambda':.02,'iterations':1000},checks={'initial_objective':float(E(pinv)),'final_objective':float(E(g)),'snr_db':snr(f,g),'relative_measurement_residual':float(np.linalg.norm(fft.fft2(g,norm='ortho')*mask-y)/np.linalg.norm(y))})
