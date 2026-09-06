from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(height=4,width=9);ax=a[0,0];diagram_axis(ax);nodes={'x':(.07,.7,r'$x$'),'y':(.34,.2,r'$y$'),'a':(.34,.7,r'$a=\log x$'),'b':(.6,.7,r'$b=ya$'),'c':(.86,.7,r'$c=\sqrt{b}$'),'f':(.86,.2,r'$f=b+c$')}
for key,(x,y,label) in nodes.items():box(ax,x,y,label,width=.16 if key=='x' else .2)
for i,j in [('x','a'),('a','b'),('y','b'),('b','c'),('c','f'),('b','f')]:
 x,y,_=nodes[i];u,v,_=nodes[j];start=(x+(.11 if u>x else 0),y-(.07 if v<y else 0));end=(u-(.11 if u>x else 0),v+(.07 if v<y else -.07 if v>y else 0));arrow(ax,start,end)
ax.text(.1,.04,r'$x>0,\quad y\log x>0$',fontsize=11)
x=2.;y=3.;b=y*np.log(x);factor=1+1/(2*np.sqrt(b));dx=y/x*factor;dy=np.log(x)*factor
finish(fig,__file__,'Exactly represents f(x,y)=y log(x)+sqrt(y log(x)) with the common product b reused by the square-root and final-addition nodes. The direct b-to-f edge is retained, which is essential for correct reverse accumulation. The differentiable domain is stated.',checks={'example_x':x,'example_y':y,'analytic_df_dx':float(dx),'analytic_df_dy':float(dy)})
