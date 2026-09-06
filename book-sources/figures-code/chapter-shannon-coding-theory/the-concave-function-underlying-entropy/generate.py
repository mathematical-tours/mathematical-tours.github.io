from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


u=np.linspace(0,1,1201);h=-special.xlogy(u,u)/np.log(2);a=1/np.e;b=1/(np.e*np.log(2))
fig,axs=canvas(height=3.5,width=6);ax=axs[0,0];curve(ax,u,h,'Entropy contribution',r'$u$',r'$h(u)=-u\log_2u$');ax.plot(a,b,'o',color=TEAL);ax.plot([a,a],[0,b],'--',color=MUTED,lw=.8);ax.annotate(r'$u=e^{-1},\quad h(u)=(e\ln2)^{-1}$',(a,b),(.5,.43),arrowprops={'arrowstyle':'->','color':MUTED});ax.set(xlim=(0,1),ylim=(0,.61));ax.scatter([0,1],[0,0],s=20,color=TEAL)
finish(fig,__file__,'Plots h(u)=-u log_2 u with its continuous value zero at u=0. The interior maximum is marked at u=1/e and height 1/(e ln 2); the curve is strictly concave for u>0.',parameters={'maximum_location':a,'maximum_value':b},checks={'nonnegative':bool(h.min()>=0)})
