from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
fig,axs=canvas(height=2.9,width=10);ax=axs[0,0];diagram_axis(ax)
xs=np.linspace(.08,.92,4)
for k,x in enumerate(xs):
    box(ax,x,.76,rf'$a_{{J+{k}}}$' if k else r'$a_J$',width=.14,height=.19)
    if k:
        box(ax,x,.25,rf'$d_{{J+{k}}}$',width=.14,height=.19)
        arrow(ax,(xs[k-1]+.08,.76),(x-.08,.76));arrow(ax,(xs[k-1]+.07,.67),(x-.07,.35))
ax.text(.5,.01,r'Approximation row $a_j$; extracted details $d_j$ below',ha='center',fontsize=11)
finish(fig,__file__,'Approximation outputs a_j are aligned on the top row; each extracted d_j is directly below the approximation from the same split. Only the approximation is recursively decomposed, and every detail band is stored once.')
