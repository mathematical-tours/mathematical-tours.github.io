from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


from matplotlib.patches import Ellipse
fig,a=canvas(1,2,height=3.8);t,f=signal_1d();a[0,0].plot(t,f)
for q in [.28,.62,.83]:a[0,0].axvspan(q-.025,q+.025,color=RED,alpha=.2)
a[0,0].set(title='Jump neighborhoods and smooth intervals',xlabel='t');f=cartoon();image(a[0,1],f,'Smooth regions and an edge curve');a[0,1].add_patch(Ellipse((.48*256,.51*256),.56*256,.68*256,fill=False,color=RED,lw=2));a[0,1].text(125,132,'regular',color='white',ha='center')
finish(fig,__file__,'Red marks identify neighborhoods of singular points in one dimension and a discontinuity curve in two dimensions. The underlying functions are piecewise smooth; the shaded jump neighborhoods represent supports that intersect singularities, not an additional discontinuity of the function.',parameters={'signal_jumps':[.28,.62,.83],'size':256})
