from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from pde import heat_example
from matplotlib.colors import to_hex
times=[0,.0001,.001,.01,.1];x,solutions=heat_example(512,times)
fig,ax=plt.subplots(figsize=(8,3.4))
for i,(time,y) in enumerate(zip(times,solutions)):
    ax.plot(np.r_[x,1],np.r_[y,y[0]],color=plt.cm.viridis(.08+.8*i/4),lw=1.1,label=f'$t={time:g}$')
ax.axhline(solutions[0].mean(),color=RED,lw=.7,ls=':',label='mean')
axes_box(ax,'$x$','$u(t,x)$');ax.set(xlim=(0,1),ylim=(-.045,1.06));ax.legend(frameon=False,fontsize=8,ncol=2,loc='upper right')
save(fig,'chapter-fourier-applications','heat-evolution')
