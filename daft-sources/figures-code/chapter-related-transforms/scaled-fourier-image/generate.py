from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from chirp import scaled_fourier2d
n=128;x=np.arange(n)-n/2;X,Y=np.meshgrid(x,x);f=((X*X+Y*Y)<16**2).astype(float)
fig,axes=plt.subplots(3,3,figsize=(7.5,7.5));axes=axes.ravel()
axes[0].imshow(1-f,cmap='gray',vmin=0,vmax=1);axes[0].set_title('Input disk')
for ax,alpha in zip(axes[1:],np.arange(.6,2.01,.2)):
    F=scaled_fourier2d(f,alpha);contrast=np.log1p(100*abs(F)/f.sum())/np.log(101)
    ax.imshow(1-contrast,cmap='gray',vmin=0,vmax=1);ax.set_title(rf'$\alpha={alpha:.1f}$')
for ax in axes:axes_box(ax);ax.set(xticks=[],yticks=[])
save(fig,'chapter-related-transforms','scaled-fourier-image')
