from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


from curvelet_frame import windows
bank,info=windows(256);idx=next(k for k,d in enumerate(info) if d['radius']==32 and abs(d['angle']-np.pi/4)<1e-8);w=bank[idx];atom=fft.fftshift(fft.ifft2(w,norm='ortho').real);spectrum=fft.fftshift(abs(fft.fft2(fft.ifftshift(atom),norm='ortho')));fig,a=canvas(1,2,height=4)
image(a[0,0],atom[96:160,96:160],'Real curvelet: spatial localization',signed=True);image(a[0,1],spectrum,'Two conjugate frequency wedges',vmin=0,vmax=spectrum.max(),cmap='viridis')
finish(fig,__file__,'The real atom is generated from a smooth compact annular wedge with parabolic radial/transverse bandwidths, and its Fourier image is computed from that same atom. Taking the real part introduces the conjugate wedge. It belongs to the explicit periodic oversampled curvelet construction used for Figure 5.26; it is not an arbitrary Gaussian labeled a tight-frame curvelet.',parameters={'size':256,'radius':32,'frequency_normal_angle':float(np.pi/4),'spatial_crop':[96,160]})
