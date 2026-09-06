from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


t,f=signal_1d(1024);fig,axs=canvas(2,2,height=5.6)
errors=[]
for ax,bins in zip(axs.ravel(),[64,32,16,8]):
    p=np.repeat(f.reshape(bins,-1).mean(axis=1),len(f)//bins);errors.append(float(np.mean((f-p)**2)))
    ax.plot(t,f,color=MUTED,lw=.8,label='signal');ax.step(t,p,where='post',label='projection');ax.set(title=rf'$j={-int(np.log2(bins))}$: {bins} intervals',xlabel=r'$t$');ax.legend(loc='lower left')
assert np.all(np.diff(errors)>=-1e-12)
finish(fig,__file__,'Each Haar projection is the exact average on its dyadic interval. Increasing j coarsens the approximation: the number of intervals decreases from 64 to 8. The same fine signal is shown in gray throughout, and the squared projection error is checked to increase when intervals are merged.',parameters={'samples':1024,'intervals':[64,32,16,8]},checks={'squared_errors':errors})
