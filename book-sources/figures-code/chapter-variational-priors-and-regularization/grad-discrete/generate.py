from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from variational_models import *


f=flower(128);d=gradient(f);fig,a=canvas(1,3,height=3.6)
for ax,(lo,hi,step) in zip(a.ravel(),[(0,128,8),(24,72,3),(60,100,2)]):
 image(ax,f[lo:hi,lo:hi],f'Pixels [{lo}, {hi})');yy,xx=np.mgrid[lo:hi:step,lo:hi:step];ax.quiver(xx-lo,yy-lo,d[1,yy,xx],d[0,yy,xx],color=RED,angles='xy',scale_units='xy',scale=.045,width=.007)
finish(fig,__file__,'Vectors use the actual forward finite differences, with columns drawn horizontally and rows vertically in image coordinates. The arrows point toward increasing intensity; the two detail panels recompute no new field and retain the same physical arrow scale.',data_sources=['data/flower.png'],parameters={'size':128,'boundary':'periodic','quiver_scale':.045})
