from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(192);rng=np.random.default_rng(20260906);versions=[f+.08*rng.standard_normal(f.shape),rng.poisson(30*f)/30,f*rng.gamma(2,.5,f.shape)]
fig,axs=canvas(1,3,height=3.4)
for ax,g,name in zip(axs.ravel(),versions,['Gaussian sensor model','Poisson photon-count model','Gamma speckle model']):image(ax,g,name)
finish(fig,__file__,'These are controlled simulations of three acquisition-noise models on the same Flower photograph, not newly measured camera, microscopy, or radar data. The Gaussian standard deviation is 0.08, Poisson means are 30 times the normalized intensity, and the mean-one Gamma multiplier has shape 2. Display clipping does not modify the simulated arrays.',data_sources=['data/flower.png'],parameters={'seed':20260906,'gaussian_sigma':.08,'poisson_peak':30,'gamma_looks':2})
