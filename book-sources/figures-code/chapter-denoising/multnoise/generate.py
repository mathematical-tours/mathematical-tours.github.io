from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=.05+.90*flower(192);rng=np.random.default_rng(20260906);looks=[1,2,4,16];fig,axs=canvas(1,4,height=3.5)
for ax,K in zip(axs.ravel(),looks):image(ax,f*rng.gamma(K,1/K,f.shape),f'K = {K} looks\nσ = {1/np.sqrt(K):.2f}')
finish(fig,__file__,'Multipliers have the explicit Gamma(shape=K,scale=1/K) law, hence mean one and variance 1/K. The Flower intensity is kept strictly positive for compatibility with the subsequent log transform. All panels use one display range; saturation is a display effect rather than truncation of the simulated model.',data_sources=['data/flower.png'],parameters={'seed':20260906,'looks':looks,'clean_intensity':'.05 + .90 flower_luminance'})
