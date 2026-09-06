from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


f=flower(192);f/=f.max();rng=np.random.default_rng(20260906);peaks=[2,8,32,128];fig,axs=canvas(1,4,height=3.4)
for ax,peak in zip(axs.ravel(),peaks):y=rng.poisson(peak*f)/peak;image(ax,y,rf'$\lambda_{{\max}}={peak}$')
finish(fig,__file__,'Replaces unspecified legacy count levels with a fully specified experiment. Photon counts are independent Poisson variables with mean lambda_max times the normalized Flower image; displayed counts are divided by lambda_max. Values above one are clipped only for display. The four acquisition parameters are now known and reproducible.',data_sources=['data/flower.png'],parameters={'seed':20260906,'maximum_mean_counts':peaks,'normalization':'Flower image divided by its maximum'})
