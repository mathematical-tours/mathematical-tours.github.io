from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


from scipy.spatial import Voronoi,voronoi_plot_2d
n=90;yy,xx=np.mgrid[:n,:n]/(n-1);points=np.c_[xx.ravel(),yy.ravel()];rng=np.random.default_rng(20260906);start=rng.random((28,2));densities=[np.ones((n,n)),.15+np.exp(-((xx-.3)**2+(yy-.65)**2)/.035)+.8*np.exp(-((xx-.75)**2+(yy-.25)**2)/.05)];fig,a=canvas(2,5,height=5.4,width=13)
for row,rho in enumerate(densities):
 c=start.copy();weights=rho.ravel().copy();weights/=weights.sum()
 for k in range(31):
  dist=np.sum((points[:,None]-c[None,:])**2,axis=2);label=dist.argmin(axis=1)
  if k in [0,2,3,5,30]:
   ax=a[row,[0,2,3,5,30].index(k)];ax.imshow(rho,extent=[0,1,0,1],origin='lower',cmap='Greys',vmin=0,vmax=2,alpha=.5);voronoi_plot_2d(Voronoi(c),ax=ax,show_points=False,show_vertices=False,line_colors=MUTED,line_width=.8,point_size=0);ax.scatter(*c.T,s=12,color=RED);ax.set(aspect='equal',title=f'{"Uniform" if row==0 else "Nonuniform"}: iteration {k}',xticks=[],yticks=[],xlim=(0,1),ylim=(0,1))
  c=np.array([np.average(points[label==j],axis=0,weights=weights[label==j]) if np.any(label==j) else c[j] for j in range(len(c))])
finish(fig,__file__,'Continuous Lloyd updates are approximated by weighted quadrature on a 90 by 90 grid. Both rows start from the same 28 centers. Nonuniform centroids use density weights, not the unweighted mean of the Voronoi pixels. Cell boundaries are the exact Euclidean Voronoi edges of the current centroids.',parameters={'seed':20260906,'quadrature_side':90,'centers':28,'display_iterations':[0,2,3,5,30]})
