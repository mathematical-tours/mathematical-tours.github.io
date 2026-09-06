from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(height=3.8);ax=axs[0,0];diagram_axis(ax)
for y,left,filt in [(.72,r'$a_{j+1}$',r'$\star h$'),(.28,r'$d_{j+1}$',r'$\star g$')]:
    box(ax,.09,y,left);box(ax,.31,y,r'$\uparrow2$');box(ax,.53,y,filt);arrow(ax,(.17,y),(.23,y));arrow(ax,(.39,y),(.45,y));arrow(ax,(.61,y),(.76,.5))
box(ax,.80,.50,'+',width=.08);box(ax,.95,.5,r'$a_j$',width=.09);arrow(ax,(.85,.5),(.90,.5))
ax.text(.5,.04,r'$a_j=(a_{j+1}\uparrow2)\star h+(d_{j+1}\uparrow2)\star g$',ha='center')
finish(fig,__file__,'The synthesis stage inserts zeros before convolution with h and g and sums the branches. These filters are not reversed in synthesis. The formula is the exact adjoint of the forward filtering/downsampling step for the real orthogonal filter convention used in the chapter.')
