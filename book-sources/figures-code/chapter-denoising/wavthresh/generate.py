from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

f,y=denoising_case();a,s=wavelet_array(y,'db4',4);T=2.7*.08
b=hard(a,T);b[s[0]]=a[s[0]];g=wavelet_inverse(b,s,'db4')
fig,axs=canvas(1,4,height=3.7,width=12)
image(axs[0,0],y,'Noisy image $f$')
wavelet_image(axs[0,1],a,s,'Wavelet coefficients')
wavelet_image(axs[0,2],b,s,'Hard thresholding')
image(axs[0,3],g,f'Reconstruction: {snr(f,g):.1f} dB')
assert np.all(b[s[1]['dd']][abs(a[s[1]['dd']])<=T]==0)
finish(fig,__file__,'Restores the original image-to-coefficients-to-thresholded-coefficients-to-image sequence using the supplied flower photograph. All detail bands use exact hard thresholding, with equality mapped to zero; the coarse band is retained. Both coefficient displays use identical scaling and boundaries.',parameters={'seed':20260906,'sigma':.08,'threshold_sigma':2.7,'wavelet':'db4','levels':4},checks={'snr_db':snr(f,g)})
