from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
t,f=signal_1d();z=np.random.default_rng(20260906).standard_normal(len(f));fig,axs=canvas(1,3,height=3.1,width=11);arrays={'clean':f}
for ax,sigma,name in zip(axs.ravel(),[.04,.12,.24],['Low','Medium','High']):
    y=f+sigma*z;ax.plot(t,f,color=MUTED,lw=.8,label=r'$f_0$');ax.plot(t,y,lw=.7,label=r'$f=f_0+w$');ax.set(title=rf'{name} noise: $\sigma={sigma:g}$',xlabel=r'$t$',ylim=(-.6,1.7));ax.legend(fontsize=8);arrays[str(sigma)]=y
finish(fig,__file__,'Three noise levels use the same clean signal and the same standard Gaussian realization, scaled by sigma. This isolates noise amplitude from changes in the underlying sample. All panels share their vertical range.',parameters={'samples':1024,'sigma':[.04,.12,.24],'seed':20260906},arrays=arrays)
