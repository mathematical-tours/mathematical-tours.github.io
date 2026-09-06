from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(1,4,height=3.3,width=13)
p={0:.27,1:.53,2:.19,3:.01};codes=huffman_codes(p)
steps=[([.01,.19,.27,.53],['$s_3$','$s_2$','$s_0$','$s_1$']),([.20,.27,.53],['$s_3,s_2$','$s_0$','$s_1$']),([.47,.53],['$s_3,s_2,s_0$','$s_1$'])]
for ax,(weights,labels),caption in zip(axs.ravel()[:3],steps,['1. Merge 0.01 + 0.19','2. Merge 0.20 + 0.27','3. Merge 0.47 + 0.53']):
    bars=ax.bar(range(len(weights)),weights,color=[RED if i<2 else TEAL for i in range(len(weights))],width=.6)
    ax.set(title=caption,ylim=(0,.66),xticks=range(len(weights)),xticklabels=labels,ylabel='Probability mass')
    for bar,value in zip(bars,weights):ax.text(bar.get_x()+bar.get_width()/2,value+.015,f'{value:.2f}',ha='center',fontsize=10)
code_tree(axs[0,3],codes,'Resulting optimal prefix tree')
H=-sum(v*np.log2(v) for v in p.values());L=sum(p[s]*len(c) for s,c in codes.items())
finish(fig,__file__,'The merges are computed from probabilities (0.27,0.53,0.19,0.01), rather than inferred from the old drawing. Red bars identify the two least likely current nodes. The resulting tree gives code lengths (2,1,3,3), expected length 1.67 bits, and the correct entropy bound.',parameters={'probabilities':p,'codes':codes},checks={'entropy':H,'expected_length':L,'entropy_bound':bool(H<=L<H+1)})
