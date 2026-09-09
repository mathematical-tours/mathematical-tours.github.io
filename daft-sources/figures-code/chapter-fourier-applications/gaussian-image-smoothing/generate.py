from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import flower
from image_processing import smooth_image
f=flower(128);k=np.arange(-24,25);X,Y=np.meshgrid(k,k)
fig,axes=plt.subplots(2,4,figsize=(9,4.5))
for j,sigma in enumerate([.3,1,2.5,6]):
    g=np.exp(-(X*X+Y*Y)/(2*sigma*sigma));g/=g.sum()
    axes[0,j].imshow(g/g.max(),cmap='gray',vmin=0,vmax=1,extent=(-24.5,24.5,24.5,-24.5))
    axes[0,j].set_title(rf'$\sigma={sigma:g}$ pixels')
    axes[1,j].imshow(smooth_image(f,sigma),cmap='gray',vmin=0,vmax=1)
    for ax in axes[:,j]:axes_box(ax);ax.set(xticks=[],yticks=[])
axes[0,0].set_ylabel('kernel / peak');axes[1,0].set_ylabel('smoothed image')
save(fig,'chapter-fourier-applications','gaussian-image-smoothing')
