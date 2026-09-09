from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from duality_data import walsh_matrix
from transforms import walsh_sequency
from matplotlib.colors import ListedColormap
from itertools import product
W=walsh_matrix(8);s=walsh_sequency(8)
order=sorted(product(range(8),repeat=2),key=lambda t:(s[t[0]]+s[t[1]],s[t[0]]))
fig,axs=plt.subplots(8,8,figsize=(6.0,6.0),layout="none")
fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.99,wspace=.13,hspace=.13)
for ax,(i,j) in zip(axs.ravel(),order):
    ax.imshow(W[i,:,None]*W[j,None,:],cmap=ListedColormap([TEAL,'white']),vmin=-1,vmax=1,interpolation='nearest')
    ax.set(xticks=[],yticks=[])
    ax.set_xticks(np.arange(.5,8),minor=True);ax.set_yticks(np.arange(.5,8),minor=True)
    ax.grid(which='minor',color=MUTED,alpha=.35,lw=.3);ax.tick_params(which='minor',length=0)
    axes_box(ax)
save(fig,'chapter-duality','walsh-basis-two-dimensions')
