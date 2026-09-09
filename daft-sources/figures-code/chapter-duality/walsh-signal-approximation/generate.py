from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from duality_data import chirp
from transforms import fwt,walsh_sequency
x,f=chirp();N=len(f);s=fwt(f);order=np.argsort(walsh_sequency(N))
fig,axs=plt.subplots(3,3,figsize=(7.1,5.0),sharex=True,sharey=True)
for ax,M in zip(axs.ravel(),[256,154,102,77,64,51,38,26,13]):
    kept=np.zeros(N,complex);kept[order[:M]]=s[order[:M]]
    approx=(fwt(kept)/N).real
    ax.plot(x,f,color=MUTED,lw=.65,alpha=.55)
    ax.plot(x,approx,color=TEAL,lw=1.05)
    ax.set(title=f'$M={M}$ ({100*M/N:.1f}%)',xlim=(0,1),ylim=(-1,1),xticks=[0,.5,1],yticks=[-1,0,1])
    axes_box(ax)
for ax in axs[-1]:ax.set_xlabel('$t$')
save(fig,'chapter-duality','walsh-signal-approximation')
