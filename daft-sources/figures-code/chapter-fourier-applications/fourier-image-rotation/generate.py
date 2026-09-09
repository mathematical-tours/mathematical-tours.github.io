from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import flower
from image_processing import rotation_examples
f=flower(97)[::-1,:].T;angles=np.linspace(0,np.pi/4,6);rotated=rotation_examples(f,angles)
center=rotated.shape[1]//2;radius=72
fig,axes=plt.subplots(1,6,figsize=(11,2.25))
for ax,y,angle in zip(axes,rotated,angles):
    cropped=y[center-radius:center+radius+1,center-radius:center+radius+1]
    ax.imshow(cropped.T,cmap='gray',vmin=0,vmax=1,origin='lower')
    axes_box(ax);ax.set(xticks=[],yticks=[],title=rf'${np.rad2deg(angle):g}^\circ$')
save(fig,'chapter-fourier-applications','fourier-image-rotation')
