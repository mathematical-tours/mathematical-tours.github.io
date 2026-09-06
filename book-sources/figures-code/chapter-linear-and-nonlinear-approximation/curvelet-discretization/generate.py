from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


from matplotlib.patches import Ellipse
fig,a=canvas(1,3,height=3.6)
for ax,(j,theta) in zip(a.ravel(),[(-4,0),(-6,0),(-6,np.pi/4)]):
 length=2**(j/2);width=2.**j;R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]);m=np.array([(k,l) for k in range(-12,13) for l in range(-40,41)]);p=(m*np.array([length,width]))@R.T+.5;ok=np.all((p>=0)&(p<=1),axis=1);ax.scatter(p[ok,0],p[ok,1],s=6,color=TEAL);ax.add_patch(Ellipse((.5,.5),length,width,angle=theta*180/np.pi,fill=False,edgecolor=RED,lw=2));ax.set(aspect='equal',xlim=(0,1),ylim=(0,1),title=rf'$j={j},\ \theta={theta/np.pi:g}\pi$');ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
finish(fig,__file__,'Centers follow exactly u=R_theta(2^(j/2)m1,2^j m2), translated to the middle of the plotting window. The highlighted support has width equal to length squared. Refining scale changes the two sampling intervals differently; rotation rotates the whole lattice.',parameters={'scales':[-4,-6,-6],'angles':[0,0,float(np.pi/4)]})
