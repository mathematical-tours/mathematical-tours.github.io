from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


fig,a=canvas(height=4,width=8);ax=a[0,0];diagram_axis(ax);layers=[]
for x,n in zip([.1,.36,.63,.9],[3,5,4,2]):layers.append([(x,y) for y in np.linspace(.2,.82,n)])
for left,right in zip(layers[:-1],layers[1:]):
 for u in left:
  for v in right:ax.plot([u[0],v[0]],[u[1],v[1]],color=MUTED,lw=.6,alpha=.4,zorder=1)
for layer in layers:
 for x,y in layer:ax.scatter(x,y,s=180,color='white',edgecolor=TEAL,zorder=3)
for k,x in enumerate([.23,.495,.765]):ax.text(x,.96,rf'$\theta_{k}$',ha='center',color=RED)
ax.text(.5,.035,r'$x_k=\rho(\theta_{k-1}x_{k-1})$ (coordinatewise activation)',ha='center',fontsize=12)
finish(fig,__file__,'A fully connected multilayer perceptron without biases shows distinct weight matrices for successive layers. Edges represent learned scalar weights, and each noninput node applies the same coordinatewise activation after summation. No weight sharing across layers is implied.',parameters={'layer_dimensions':[3,5,4,2],'biases':False})
