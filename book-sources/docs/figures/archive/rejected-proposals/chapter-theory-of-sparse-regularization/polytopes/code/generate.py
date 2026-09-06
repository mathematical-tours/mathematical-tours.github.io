from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.optimize import linprog
A=np.array([[1,.2,-.7],[0,.9,.4]]);vertices=np.r_[np.eye(3),-np.eye(3)];Y=vertices@A.T;hull=ConvexHull(Y);fig=plt.figure(figsize=(9,4),layout='constrained');ax=fig.add_subplot(121,projection='3d');faces=ConvexHull(vertices).simplices;ax.add_collection3d(Poly3DCollection(vertices[faces],facecolor=TEAL,alpha=.15,edgecolor=TEAL));ax.set(xlim=(-1,1),ylim=(-1,1),zlim=(-1,1),title=r'$B_1\subset\mathbb{R}^3$',xlabel='x₁',ylabel='x₂',zlabel='x₃');ax=fig.add_subplot(122);poly=Y[hull.vertices];ax.fill(poly[:,0],poly[:,1],color=TEAL,alpha=.15);ax.plot(*np.vstack([poly,poly[0]]).T,color=TEAL);ax.scatter(*Y.T,color=TEAL);x0=np.array([.6,0,0]);y=A@x0;sol=linprog(np.ones(6),A_eq=np.c_[A,-A],b_eq=y,bounds=(0,None),method='highs');assert sol.success;recovered=sol.x[:3]-sol.x[3:];ax.scatter(*y,color=RED);ax.set(aspect='equal',title=r'$A(B_1)\subset\mathbb{R}^2$',xlabel='y₁',ylabel='y₂');fig.suptitle('Linear projection A; recovery chooses a minimum-ℓ¹ preimage',fontsize=12)
finish(fig,__file__,'The octahedral ell-1 ball is mapped by an explicit 2 by 3 matrix; the right polygon is its computed convex hull. A marked measurement has a minimum-ell-1 preimage computed by linear programming. The recovery map is not depicted as an inverse on all of R³.',parameters={'A':A.tolist()},checks={'selected_preimage':recovered.tolist(),'measurement_residual':float(np.linalg.norm(A@recovered-y))})
