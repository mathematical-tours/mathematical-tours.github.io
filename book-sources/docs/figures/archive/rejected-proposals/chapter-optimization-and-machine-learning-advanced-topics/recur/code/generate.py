from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(height=3.8,width=9);ax=a[0,0];diagram_axis(ax);xs=[.07,.29,.51,.73]
for k,x in enumerate(xs):box(ax,x,.34,rf'$x_{k}$',width=.12)
for x in xs[1:]:arrow(ax,(x-.15,.34),(x-.07,.34),r'$h$')
box(ax,.94,.34,r'$L(x_3,\theta)$',width=.16,size=11);arrow(ax,(.80,.34),(.85,.34));box(ax,.5,.84,r'$\theta$',width=.16,color=RED)
for x in xs[1:]+[.94]:arrow(ax,(.5,.75),(x,.44),color=RED)
ax.text(.5,.04,'One shared parameter block; include its direct effect on the loss',ha='center',fontsize=11)
finish(fig,__file__,'One parameter node feeds every recurrence and the terminal loss, matching f(theta)=L(x_t,theta). Its gradient sums contributions from all unrolled steps plus the explicit loss derivative. Distinct copies of theta are not drawn as independent parameters.')
