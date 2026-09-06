from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


f,H,y,grid,runs=deblur_sweep();fig,a=canvas(2,3,height=6);image(a[0,0],f,'Original');image(a[0,1],y,f'Observation: {snr(f,y):.1f} dB');records={}
for ax,(name,images) in zip([a[0,2],a[1,0],a[1,1]],runs.items()):
 scores=[snr(f,g) for g in images];k=int(np.argmax(scores));image(ax,images[k],f'{name}\n{scores[k]:.1f} dB');records[name]={'lambda':float(grid[k]),'snr_db':scores[k]}
lam=records['Orthogonal wavelets']['lambda'];g=sparse_deblur(y,H,lam,True,600);image(a[1,2],g,f'Stationary synthesis frame\n{snr(f,g):.1f} dB');W,S=wavelet_operator(f.shape,True);rng=np.random.default_rng(9);c=rng.standard_normal(W(f).shape);err=abs(np.sum(W(f)*c)-np.sum(f*S(c)));assert err<1e-9
finish(fig,__file__,'Squared norm, Sobolev and orthogonal wavelet reconstructions use the exact grid plotted in Figure 10.9. The last panel runs coefficient-space synthesis FISTA with a normalized stationary Parseval frame; it does not pretend that analysis-threshold-synthesis is the proximal map of an analysis penalty. Each stationary coefficient is penalized by its atom L2 norm, equivalently using a unit-norm synthesis dictionary. Its lambda equals the orthogonal oracle choice and is stated as such. The frame adjoint identity is checked.',data_sources=['data/flower.png'],parameters={'size':256,'seed':20260906,'blur_sigma':1.2,'noise_sigma':.015,'stationary_lambda':lam},checks={'oracle_choices':records,'stationary_adjoint_error':float(err),'stationary_snr_db':snr(f,g)})
