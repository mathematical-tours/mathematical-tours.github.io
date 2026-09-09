from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from duality_data import walsh_matrix
from transforms import walsh_sequency
from matplotlib.colors import ListedColormap
W=walsh_matrix(64);order=np.argsort(walsh_sequency(64))
fig,axs=plt.subplots(1,2,figsize=(7.1,3.35))
for ax,a,title,ylabel in zip(axs,[W,W[order]],['Natural order','Sequency order'],['Character index','Sign changes']):
    ax.imshow(a,cmap=ListedColormap([TEAL,'white']),vmin=-1,vmax=1,interpolation='nearest')
    ax.set(title=title,xticks=[0,16,32,48,63],yticks=[0,16,32,48,63])
    axes_box(ax,'Sample index',ylabel)
save(fig,'chapter-duality','walsh-orderings')
