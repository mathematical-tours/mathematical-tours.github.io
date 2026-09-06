from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *

f=flower(256);row,col=np.indices(f.shape)
mask=((col%40)>=5)&((row%64)>=4);y=f*mask
g=inpaint(y,mask,'sobolev',1500);fig,axs=canvas(1,3,height=3.8)
for ax,z,title in zip(axs.ravel(),[f,y,g],['Supplied flower photograph','Cage-shaped missing-pixel mask','Sobolev inpainting']):image(ax,z,title)
assert np.max(abs(g[mask]-f[mask]))<1e-12
finish(fig,__file__,'Recreates the cage-removal setup using the supplied flower photograph and an explicit grid-shaped occlusion. Both the known-pixel mask and the clean reference are available for this controlled experiment. The Sobolev reconstruction enforces observed pixels exactly after every step; no legacy parrot image is used.',parameters={'size':256,'vertical_period':40,'vertical_width':5,'horizontal_period':64,'horizontal_width':4,'iterations':1500,'step':.24},checks={'known_pixel_error':float(np.max(abs(g[mask]-y[mask]))),'observed_fraction':float(mask.mean()),'snr_db':snr(f,g)})
