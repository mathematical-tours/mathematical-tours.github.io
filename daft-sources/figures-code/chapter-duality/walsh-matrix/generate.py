from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from duality_data import walsh_matrix
from matplotlib.colors import ListedColormap
W=walsh_matrix(64)
assert np.array_equal(W@W.T,64*np.eye(64,dtype=int))
fig,ax=plt.subplots(figsize=(4.0,3.9))
ax.imshow(W,cmap=ListedColormap([TEAL,'white']),vmin=-1,vmax=1,interpolation='nearest')
ax.set(xticks=[0,16,32,48,63],yticks=[0,16,32,48,63])
axes_box(ax,'Sample index $j$','Character index $i$')
save(fig,'chapter-duality','walsh-matrix')
