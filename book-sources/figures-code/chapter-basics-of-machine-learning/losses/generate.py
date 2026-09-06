from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


s=np.linspace(-3,2,1200);fig,a=canvas(height=4,width=7);ax=a[0,0];ax.plot(s,np.maximum(1+s,0),label='hinge');ax.plot(s,np.logaddexp(0,s),label='logistic (unscaled)');ax.plot(s,np.logaddexp(0,s)/np.log(2),'--',label='logistic / log 2');ax.plot(s,np.exp(s),label='exponential');ax.plot(s,(1+s)**2,label='squared');ax.plot([-3,0],[0,0],color=INK,lw=2);ax.plot([0,2],[1,1],color=INK,lw=2,label='0–1');ax.scatter([0],[0],facecolor='white',edgecolor=INK,zorder=5);ax.scatter([0],[1],color=INK,zorder=5);ax.set(ylim=(-.05,4),xlabel=r'$s=-y\langle x,\beta\rangle$',ylabel='Loss',title='Classification losses versus negative margin');ax.legend(ncol=2);ax.grid(True)
finish(fig,__file__,'Losses are evaluated exactly as functions of negative margin. The binary loss equals one at zero; open/closed dots show the discontinuity. Both logistic normalizations are distinguished, and division by log(2) makes the binary-loss upper bound tight at zero. Vertical truncation is a stated viewing range.',parameters={'s_range':[-3,2],'vertical_display_limit':4},checks={'normalized_logistic_at_zero':float(np.logaddexp(0,0)/np.log(2))})
