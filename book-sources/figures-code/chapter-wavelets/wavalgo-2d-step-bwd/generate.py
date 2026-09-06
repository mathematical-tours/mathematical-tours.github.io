from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(height=4.7);ax=axs[0,0];diagram_axis(ax)
for yy,y,label,g in [(.74,.88,r'$a_j$','h'),(.74,.62,r'$d_j^V$','g'),(.25,.36,r'$d_j^H$','h'),(.25,.10,r'$d_j^D$','g')]:
    box(ax,.07,y,label,width=.12);arrow(ax,(.14,y),(.44,yy),rf'$\uparrow^V2,\ \star^V {g}$')
for y,label,h in [(.74,r'$\tilde a_j$','h'),(.25,r'$\tilde d_j$','g')]:
    box(ax,.5,y,label,width=.12);ax.text(.41,y,'+',ha='center',va='center');arrow(ax,(.57,y),(.88,.5),rf'$\uparrow^H2,\ \star^H {h}$')
box(ax,.94,.5,r'$a_{j-1}$',width=.11);ax.text(.85,.5,'+',ha='center',va='center');ax.text(.5,.98,'Sum paired branches, then reconstruct the fine approximation',ha='center')
finish(fig,__file__,'The inverse graph reverses the analysis order: first reconstruct along the second coordinate from (a_j,d_j^V) and (d_j^H,d_j^D), then along the first coordinate. Each branch inserts zeros before convolution with the unreversed synthesis filter. Plus signs mark the two intermediate sums and the final sum.')
