from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f,y=denoising_case();ht=[0,.15,.6,2.4];tt=[0,.02,.06,.18];tvf=tv_flow(y,tt);fig,a=canvas(2,4,height=6)
for k in range(4):
 for row,g,t in [(0,heat(y,ht[k]),ht[k]),(1,tvf[k],tt[k])]:image(a[row,k],g,f'{"Heat" if row==0 else "TV"}: t = {t}\n{snr(f,g):.1f} dB')
finish(fig,__file__,'Flow snapshots use one fixed noisy flower. Heat and smoothed-TV times are stated separately, and SNR is recomputed before display clipping. TV time integration uses a stable step below epsilon/4.',data_sources=['data/flower.png'],parameters={'seed':20260906,'sigma':.08,'epsilon':.01,'heat_times':ht,'tv_times':tt})
