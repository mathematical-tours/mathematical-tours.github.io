from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


rng=np.random.default_rng(20260906);n=256;a0=np.zeros(n);support=rng.choice(n,12,replace=False);a0[support]=rng.choice([-1,1],12)*rng.uniform(.6,1.4,12);a=a0+.12*rng.standard_normal(n);fig,axs=canvas(1,2,height=3.5)
for ax,v,title in [(axs[0,0],a0,'Clean sparse coefficients'),(axs[0,1],a,'Noisy coefficients')]:ax.vlines(np.arange(n),0,v,color=TEAL,lw=.8);ax.set(title=title,xlabel='Coefficient index',ylim=(-1.7,1.7));ax.grid(True)
finish(fig,__file__,'A 256-dimensional coefficient vector with exactly twelve nonzeros is corrupted by independent Gaussian noise of standard deviation 0.12. Signal and noise are generated from a fixed seed; both plots have the same axes.',parameters={'seed':20260906,'dimension':n,'sparsity':12,'sigma':.12},checks={'nonzero_count':int(np.count_nonzero(a0))},arrays={'clean_coefficients':a0,'noisy_coefficients':a})
