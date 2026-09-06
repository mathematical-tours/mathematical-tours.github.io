from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f=flower(128);mask=np.random.default_rng(20260906).random(f.shape)>.7;y=f*mask;fig,a=canvas(1,3,height=3.7);records={};image(a[0,0],y,'30% known pixels')
for ax,kind in zip(a[0,1:],['sobolev','tv']):g=inpaint(y,mask,kind,1000);image(ax,g,f'{kind.upper()}: {snr(f,g):.1f} dB');records[kind]={'snr_db':snr(f,g),'known_pixel_error':float(np.max(abs(g[mask]-y[mask])))}
finish(fig,__file__,'Sobolev and isotropic nonsmooth-TV inpainting share the exact same random observation mask. Known pixels are imposed as equality constraints at every iteration. Recomputed SNR values determine the outcome rather than preserving the ranking of a different legacy dataset.',data_sources=['data/flower.png'],parameters={'size':128,'seed':20260906,'observation_probability':.3,'iterations':1000},checks=records)
