from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(1,2,height=4.8);ax=axs[0,0];ax.set_title('Four-row stripe scan')
order=[(r,c) for start in [0,4] for c in range(5) for r in range(start,start+4)]
for k,(r,c) in enumerate(order):
    ax.add_patch(plt.Rectangle((c,r),1,1,facecolor='#E9F2F3' if r<4 else '#F8EEDF',edgecolor='white'));ax.text(c+.5,r+.5,str(k+1),ha='center',va='center',fontsize=10)
for start in [0,4]:
    for c in range(5):ax.annotate('',xy=(c+.84,start+3.78),xytext=(c+.84,start+.16),arrowprops={'arrowstyle':'->','lw':.8,'color':TEAL})
ax.set(xlim=(0,5),ylim=(8,0),aspect='equal');ax.set_axis_off()
ax=axs[0,1];ax.set_title('Local 3 × 3 coding context')
for r in range(3):
 for c in range(3):
    center=(r,c)==(1,1);ax.add_patch(plt.Rectangle((c,r),1,1,facecolor=TEAL if center else '#E9F2F3',edgecolor='white',lw=2));ax.text(c+.5,r+.5,'current\ncoefficient' if center else 'known\nstate',ha='center',va='center',fontsize=10,color='white' if center else INK)
ax.set(xlim=(0,3),ylim=(3.9,-.1),aspect='equal');ax.set_axis_off();ax.text(1.5,3.4,'Significance · sign · refinement',ha='center',fontsize=10)
finish(fig,__file__,'The numbered scan visits stripes of height four, column by column within each stripe. The context panel distinguishes the coefficient being coded from the neighboring states already known to the decoder; it does not use future unrevealed bits. Exact standard-specific context tables are not substituted by this schematic.',parameters={'stripe_height':4,'demonstration_block':[8,5]},checks={'each_position_visited_once':len(set(order))==40})
