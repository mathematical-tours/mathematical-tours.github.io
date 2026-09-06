from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


from matplotlib.patches import Circle,Ellipse
fig,a=canvas(1,2,height=4)
for ax,curvelets in zip(a.ravel(),[False,True]):
 ax.add_patch(Circle((0,0),1,facecolor='#edf1f2',edgecolor=INK,lw=1.1))
 for theta in np.linspace(.1,2*np.pi+.1,15,endpoint=False):
  x,y=np.cos(theta),np.sin(theta);angle=(theta*180/np.pi+90) if curvelets else 90*round((theta*180/np.pi+90)/90)
  for offset in [-.035,0,.035]:
   ax.add_patch(Ellipse(((1+offset)*x,(1+offset)*y),.42 if curvelets else .19,.014 if curvelets else .045,angle=angle,facecolor=TEAL,edgecolor='none',alpha=.65))
 # A mismatched atom crossing the interface has weaker cancellation alignment.
 theta=.15;x,y=np.cos(theta),np.sin(theta);ax.add_patch(Ellipse((1.1*x,1.1*y),.43 if curvelets else .2,.022 if curvelets else .05,angle=theta*180/np.pi,facecolor=RED,edgecolor='none'))
 ax.text(0,-1.35,'Isotropic scale; fixed orientations' if not curvelets else 'Parabolic scale; tangent orientations',ha='center',fontsize=10);ax.set(aspect='equal',xlim=(-1.45,1.45),ylim=(-1.48,1.4),title='Directional wavelets' if not curvelets else 'Curvelets');ax.set_axis_off()
finish(fig,__file__,'Reconstructs the original geometric schematic rather than replacing it with an unrelated numerical panel. Wavelet footprints have comparable tangential/normal scales and a few fixed orientations; curvelet footprints are elongated and follow the contour tangent. Red transverse atoms illustrate an orientation mismatch. Footprints indicate localization, not exact compact supports or numerical coefficient values.')
