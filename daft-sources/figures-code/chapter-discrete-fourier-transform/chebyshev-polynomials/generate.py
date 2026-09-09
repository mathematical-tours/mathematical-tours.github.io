from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
x=np.linspace(-1,1,501)
fig,axes=plt.subplots(2,3,figsize=(9,4.3),sharex=True,sharey=True)
for k,ax in enumerate(axes.flat,1):
    ax.plot(x,np.cos(k*np.arccos(x)),color=TEAL);ax.axhline(0,color=MUTED,lw=.5)
    axes_box(ax,'$x$');ax.set(title=rf'$T_{k}$',xlim=(-1,1),ylim=(-1.12,1.12),yticks=[-1,0,1])
save(fig,'chapter-discrete-fourier-transform','chebyshev-polynomials')
