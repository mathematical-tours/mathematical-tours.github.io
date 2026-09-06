from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(height=3.3);ax=axs[0,0];diagram_axis(ax)
for x,y,label,w in [(.10,.42,r'$f_0$',.13),(.33,.42,'+',.07),(.51,.42,r'$Y=f_0+w$',.20),(.76,.42,r'$\mathcal{D}$',.12),(.95,.42,r'$\widehat f$',.09),(.33,.8,r'$w$',.10)]:box(ax,x,y,label,width=w)
for a,b in [((.18,.42),(.28,.42)),((.38,.42),(.40,.42)),((.62,.42),(.69,.42)),((.83,.42),(.90,.42)),((.33,.73),(.33,.49))]:arrow(ax,a,b)
ax.text(.1,.19,'clean signal',ha='center',fontsize=10);ax.text(.51,.19,'observation',ha='center',fontsize=10);ax.text(.76,.19,'denoiser',ha='center',fontsize=10);ax.text(.33,.95,'acquisition noise',ha='center',fontsize=10)
finish(fig,__file__,'The clean signal and additive noise meet at a sum node, producing the observation Y. Only Y enters the denoiser D; the clean target remains unavailable to the estimator. The reconstruction is labeled as an estimate, not an exact inverse of the noisy acquisition.')
