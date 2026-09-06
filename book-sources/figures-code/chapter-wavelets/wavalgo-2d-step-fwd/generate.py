from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(height=4.7);ax=axs[0,0];diagram_axis(ax)
box(ax,.06,.50,r'$a_{j-1}$',width=.11)
for y,label,h in [(.74,r'$\tilde a_j$','h'),(.25,r'$\tilde d_j$','g')]:
    box(ax,.45,y,label,width=.15);arrow(ax,(.12,.5),(.37,y),rf'$\star^H\bar {h},\ \downarrow^H2$')
for yy,y,label,g in [(.74,.88,r'$a_j$','h'),(.74,.62,r'$d_j^V$','g'),(.25,.36,r'$d_j^H$','h'),(.25,.10,r'$d_j^D$','g')]:
    box(ax,.92,y,label,width=.12);arrow(ax,(.53,yy),(.85,y),rf'$\star^V\bar {g},\ \downarrow^V2$')
ax.text(.25,.98,'First coordinate',ha='center');ax.text(.71,.98,'Second coordinate',ha='center')
finish(fig,__file__,'The 2-D analysis graph follows the chapter’s definitions exactly: H is filtering/downsampling in the first coordinate and V in the second. Low-low yields a_j, low-high yields d_j^V, high-low yields d_j^H, and high-high yields d_j^D. Every downsampling follows convolution with the reversed analysis filter.')
