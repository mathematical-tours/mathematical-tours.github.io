from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from duality_data import letter_a
from transforms import walsh_approx2d
x=letter_a(64);N=x.size
fig,axs=plt.subplots(3,3,figsize=(6.7,6.3))
errors=[]
for ax,M in zip(axs.ravel(),[4096,3277,2458,2048,1638,1229,819,410,205]):
    result=walsh_approx2d(x,M);errors.append(np.linalg.norm(result-x))
    ax.imshow(result,cmap='gray',vmin=0,vmax=1,interpolation='nearest')
    ax.set(title=f'$M={M}$ ({100*M/N:.0f}%)',xticks=[],yticks=[])
    axes_box(ax)
assert np.all(np.diff(errors)>=-1e-10)
save(fig,'chapter-duality','walsh-image-approximation')
