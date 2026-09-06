from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


from matplotlib.patches import Circle
angles=np.array([.45,1.1,1.8,2.6]);points=np.c_[np.cos(angles),np.sin(angles)];signed=np.r_[points,-points];theta=(.45+1.1)/2;n=np.array([np.cos(theta),np.sin(theta)]);height=np.cos((1.1-.45)/2);p=n/height;t=np.array([-n[1],n[0]]);fig,a=canvas(height=4.5,width=6);ax=a[0,0];ax.add_patch(Circle((0,0),1,facecolor='#f0f4f4',edgecolor=MUTED));xx=np.linspace(-1.15,1.15,2);line=n[:,None]*height+t[:,None]*xx;ax.plot(*line,color=RED);ax.scatter(*signed.T,color=TEAL,s=35);ax.scatter(*points[:2].T,color=RED,s=55);ax.annotate('',xy=p,xytext=(0,0),arrowprops={'arrowstyle':'->','color':INK});ax.text(*(p*1.07),r'$p$');ax.text(-.9,-1.22,r'$\langle p,\pm a_i\rangle\leq1$');ax.set(aspect='equal',xlim=(-1.4,1.4),ylim=(-1.35,1.4),title='Supporting line and an empty spherical cap');ax.set_axis_off();bound=float(np.max(signed@p));assert bound<=1+1e-12
finish(fig,__file__,'All signed columns lie on the unit circle. The line <p,a>=1 passes through two selected adjacent columns, and the open cap beyond it contains none of the signed columns. A direct dot-product check verifies the certificate rather than relying only on the drawing.',parameters={'column_angles':angles.tolist()},checks={'maximum_signed_pairing':bound,'support_pairings':(points[:2]@p).tolist()})
