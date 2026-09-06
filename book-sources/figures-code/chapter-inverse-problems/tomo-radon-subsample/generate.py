from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


from skimage.transform import radon
f=phantom(128);theta=np.arange(12)*180/12;sino=radon(f,theta=theta,circle=False);mask=radial_mask(128,12);fig,a=canvas(1,3,height=3.8);image(a[0,0],f,'Original tomography phantom');o=a[0,1].imshow(sino,cmap='gray',aspect='auto',extent=[0,180,sino.shape[0],0]);a[0,1].set(title='12 discrete Radon projections',xlabel='Normal angle (degrees)',ylabel='Detector position');image(a[0,2],fft.fftshift(mask),'Cartesian radial Fourier mask')
finish(fig,__file__,'The middle panel is an actual discrete Radon transform. The last panel illustrates the separate Cartesian radial Fourier model used for reconstruction in Figures 9.9–9.10. The Fourier slice theorem motivates this model, but rasterized radial sampling is not asserted to equal the interpolation-based discrete Radon transform exactly.',parameters={'size':128,'angles':theta.tolist(),'fourier_mask_line_halfwidth':.5})
