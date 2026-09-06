from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f=flower(128);rng=np.random.default_rng(20260906);mask=rng.random(f.shape)>.7;fig,a=canvas(1,4,height=3.5)
for ax,g,title in zip(a.ravel(),[f,ndimage.gaussian_filter(f,2,mode='wrap'),f[::4,::4],np.where(mask,f,0)],['Unknown image','Blur','Subsampling × 4','30% observed pixels']):image(ax,g,title)
finish(fig,__file__,'Forward operators are applied to a common acquired image: periodic Gaussian convolution, restriction to one pixel in each 4 by 4 cell, and a known random observation mask. These are measurements, not reconstructions.',data_sources=['data/flower.png'],parameters={'size':128,'seed':20260906,'blur_sigma':2,'subsampling':4,'observation_probability':.3})
