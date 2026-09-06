from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


f,H,y,grid,runs=deblur_sweep();fig,a=canvas(1,3,height=3.6);arrays={'lambda':grid}
for ax,(name,images) in zip(a.ravel(),runs.items()):scores=np.array([snr(f,g) for g in images]);curve(ax,grid,scores,name,r'$\lambda$','SNR (dB)');ax.set_xscale('log');k=int(np.argmax(scores));ax.scatter(grid[k],scores[k],color=RED);arrays[name]=scores
finish(fig,__file__,'Measured deconvolution SNR curves use the same observation, parameters, and solved objectives as Figure 10.8. Red markers are oracle maxima on the finite lambda grid, not an automatic parameter-selection method available without a clean reference.',data_sources=['data/flower.png'],parameters={'size':256,'seed':20260906,'noise_sigma':.015,'blur_sigma':1.2,'lambda_grid':grid.tolist()},arrays=arrays)
