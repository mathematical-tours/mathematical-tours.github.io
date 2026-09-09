from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import flower
from image_processing import smooth_image
f=flower(192);y=f+.12*np.random.default_rng(20260908).standard_normal(f.shape);z=smooth_image(y,.9)
fig,axes=plt.subplots(1,2,figsize=(7.4,3.7))
for ax,image,title in zip(axes,[y,z],['Noisy image','Gaussian smoothing']):
    ax.imshow(image,cmap='gray',vmin=0,vmax=1);axes_box(ax);ax.set(xticks=[],yticks=[],title=title)
    ax.set_xlabel(f'PSNR = {-10*np.log10(np.mean((image-f)**2)):.1f} dB')
assert np.mean((z-f)**2)<np.mean((y-f)**2)
save(fig,'chapter-fourier-applications','image-denoising')
