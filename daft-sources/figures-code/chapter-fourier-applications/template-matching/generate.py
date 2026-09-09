from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import flower
from image_processing import normalized_correlation
from scipy.signal import fftconvolve
from matplotlib.patches import Rectangle
f=flower(128);u,v,p=26,53,25;g=f[u:u+p,v:v+p]
raw=fftconvolve(f,g[::-1,::-1],mode='valid');normalized=normalized_correlation(f,g)
assert np.unravel_index(np.nanargmax(normalized),normalized.shape)==(u,v)
fig,axes=plt.subplots(2,2,figsize=(8,7.3))
axes[0,0].imshow(f,cmap='gray',vmin=0,vmax=1);axes[0,0].add_patch(Rectangle((v-.5,u-.5),p,p,fill=False,edgecolor=RED,lw=1.8));axes[0,0].set_title('Image and selected patch')
axes[0,1].imshow(g,cmap='gray',vmin=0,vmax=1);axes[0,1].set_title('Template')
im=axes[1,0].imshow(raw,cmap='magma');fig.colorbar(im,ax=axes[1,0],shrink=.78,label='raw correlation')
im=axes[1,1].imshow(normalized,cmap='RdBu_r',vmin=-1,vmax=1);fig.colorbar(im,ax=axes[1,1],shrink=.78,label='normalized correlation')
axes[1,0].set_title('Raw correlation');axes[1,1].set_title('Normalized correlation')
for ax in axes[0]:axes_box(ax);ax.set(xticks=[],yticks=[])
for ax in axes[1]:ax.plot(v,u,'x',color='#40ffbc',ms=8,mew=1.8);axes_box(ax,'patch column','patch row')
save(fig,'chapter-fourier-applications','template-matching')
