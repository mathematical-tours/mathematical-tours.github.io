from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


probabilities=[[0,1,0,0],[.25]*4,[2/25,9/25,12/25,2/25]];fig,axs=canvas(1,3,height=3.2)
entropies=[]
for ax,p in zip(axs.ravel(),probabilities):
    p=np.array(p);H=-np.sum(special.xlogy(p,p))/np.log(2);entropies.append(float(H));ax.bar(np.arange(4),p,color=TEAL,width=.65);ax.set(title=rf'$H(p)={H:.3g}$ bits',ylim=(0,1.13),xticks=np.arange(4),xticklabels=[r'$p_0$',r'$p_1$',r'$p_2$',r'$p_3$'],ylabel='Probability');
    for k,value in enumerate(p):ax.text(k,value+.035,f'{value:g}',ha='center',fontsize=10)
finish(fig,__file__,'Recomputed the three normalized distributions and their base-two entropies. The uniform four-symbol distribution has exactly two bits of entropy. All bars share a vertical probability scale; the third entropy is rounded from its computed value.',parameters={'probabilities':probabilities},checks={'entropies_bits':entropies,'normalization':all(abs(sum(p)-1)<1e-12 for p in probabilities)})
