from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(height=3.8,width=8);ax=a[0,0];diagram_axis(ax);nodes={1:(.08,.75),2:(.08,.25),3:(.35,.75),4:(.35,.25),5:(.64,.5),6:(.9,.5)}
for i,(x,y) in nodes.items():box(ax,x,y,rf'$x_{i}$',width=.12)
edges=[(1,3),(1,4),(2,4),(3,5),(4,5),(3,6),(5,6)]
for i,j in edges:
 x,y=nodes[i];u,v=nodes[j]
 if (i,j)==(3,6):ax.annotate('',xy=(u,v+.08),xytext=(x+.05,y+.07),arrowprops={'arrowstyle':'->','color':MUTED,'connectionstyle':'arc3,rad=-.25'})
 else:arrow(ax,(x+.07,y),(u-.07,v))
ax.text(.08,.04,'Inputs',ha='center');ax.text(.5,.04,'Intermediate values',ha='center');ax.text(.9,.04,'Output',ha='center');assert all(i<j for i,j in edges)
finish(fig,__file__,'A directed acyclic graph in topological order illustrates shared intermediate results and branching dependencies. Every edge points from a smaller index to a larger one; no edge implies dependence on an unrelated input.',parameters={'edges':edges},checks={'topological_order_valid':True})
