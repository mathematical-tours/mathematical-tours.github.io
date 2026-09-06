from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f,y,runs=flow_experiment();fig,a=canvas(2,3,height=6);image(a[0,0],y,f'Noisy: {snr(f,y):.1f} dB');records={}
for ax,(name,(grid,images)) in zip([a[0,1],a[0,2],a[1,0],a[1,1]],runs.items()):
 scores=[snr(f,g) for g in images];k=int(np.argmax(scores));image(ax,images[k],f'{name}\n{scores[k]:.1f} dB');records[name]={'parameter':float(grid[k]),'snr_db':scores[k]}
best,_,_=tune_threshold(f,y,.08,'hard',True);image(a[1,2],best[1],f'Stationary wavelets\n{best[0]:.1f} dB')
finish(fig,__file__,'Each method is computed from the same observation and uses the clean reference only for explicitly oracle-based parameter selection. The four variational grid searches are exactly those plotted in Figure 8.9. TV regularization solves the nonsmooth ROF objective by a primal-dual algorithm; the TV flow uses epsilon=.01.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08},checks=records)
