from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f,y=denoising_case();widths=np.linspace(.1,6,60);scores=[snr(f,ndimage.gaussian_filter(y,s,mode='wrap')) for s in widths];s=widths[np.argmax(scores)];g=ndimage.gaussian_filter(y,s,mode='wrap');fig,axs=canvas(1,2,height=4)
image(axs[0,0],y,f'Noisy: {snr(f,y):.1f} dB');image(axs[0,1],g,f'Filtered: s = {s:.2f}, {snr(f,g):.1f} dB')
finish(fig,__file__,'Uses exactly the clean/noisy flower pair, width grid, and oracle SNR maximizer in the filter-selection plot. Both panels share a fixed intensity range. The image changes are the actual output of the normalized periodic Gaussian filter.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'selected_width':float(s)},checks={'snr_db':snr(f,g)})
