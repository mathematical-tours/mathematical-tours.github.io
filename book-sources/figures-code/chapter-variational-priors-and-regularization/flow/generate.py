from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f=flower(256);ht=[0,1,4,16];tt=[0,.03,.12,.4];tvf=tv_flow(f,tt);fig,a=canvas(2,4,height=5.8)
for k in range(4):image(a[0,k],heat(f,ht[k]),f'Heat: t = {ht[k]}');image(a[1,k],tvf[k],f'Smoothed TV: t = {tt[k]}')
energies=[tv(g,.01) for g in tvf];assert np.all(np.diff(energies)<=1e-8)
finish(fig,__file__,'Exact discrete heat evolution is compared with explicit smoothed-TV gradient flow (epsilon=.01, step=.002 < epsilon/4). Each row states its own times because the two energy scalings differ. The common periodic flower image and unit display range expose diffusion across edges.',parameters={'size':256,'heat_times':ht,'tv_times':tt,'epsilon':.01,'tv_step':.002},checks={'tv_energies':energies,'mean_preservation_error':float(max(abs(g.mean()-f.mean()) for g in tvf))})
