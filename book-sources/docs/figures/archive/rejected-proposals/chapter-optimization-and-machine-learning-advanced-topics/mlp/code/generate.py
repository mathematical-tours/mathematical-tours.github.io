from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(height=3.5,width=9);ax=a[0,0];diagram_axis(ax)
for k,x in enumerate([.08,.3,.52,.74]):box(ax,x,.4,rf'$x_{k}$',width=.12)
for k,x in enumerate([.3,.52,.74]):box(ax,x,.82,rf'$\theta_{k}$',width=.14,color=RED);arrow(ax,(x,.73),(x,.49));arrow(ax,(x-.15,.4),(x-.07,.4),rf'$f_{k+1}$')
box(ax,.94,.4,r'$L(x_3)$',width=.14);arrow(ax,(.81,.4),(.86,.4));ax.text(.5,.06,r'$x_k=f_k(x_{k-1},\theta_{k-1})$',ha='center',fontsize=15)
finish(fig,__file__,'Each layer has its own parameter block theta_(k-1), used only to produce x_k from x_(k-1). The scalar terminal loss has no direct parameter dependency in this feedforward example. Layer indices agree with the chapter formulas.')
