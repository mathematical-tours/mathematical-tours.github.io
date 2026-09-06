from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

f=flower(256)
a,s=wavelet_array(f,'db2',4)
fig,axs=canvas(1,2,height=4)
image(axs[0,0],f,'Image $f$')
wavelet_image(axs[0,1],a,s,'Wavelet coefficients',labels=True)
err=np.max(abs(wavelet_inverse(a,s,'db2')-f));assert err<1e-12
finish(fig,__file__,'Uses the supplied flower image while restoring the original coefficient layout. The orthonormal Daubechies-2 transform uses periodic boundaries and four decomposition levels. The coarse image and three orientations at each scale occupy their exact packed locations.',parameters={'wavelet':'db2','levels':4},checks={'inverse_max_error':float(err)})
