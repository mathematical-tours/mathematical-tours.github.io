from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


from sklearn.neighbors import KNeighborsClassifier
from matplotlib.colors import ListedColormap
X,y=iris();m,s,V,Z=pca(X);Z=Z[:,:2];xx,yy,grid=grid2(Z,200,.7);fig,a=canvas(1,3,height=3.8)
for ax,R in zip(a.ravel(),[1,5,25]):
 model=KNeighborsClassifier(R).fit(Z,y);pred=model.predict(grid).reshape(xx.shape);ax.imshow(pred,origin='lower',extent=[xx.min(),xx.max(),yy.min(),yy.max()],cmap=ListedColormap(COLORS[:3]),alpha=.15,vmin=0,vmax=2,aspect='auto',interpolation='nearest')
 for k in range(3):ax.scatter(*Z[y==k].T,color=COLORS[k],s=10)
 ax.set(title=f'R = {R} neighbors',xlabel='PC1',ylabel='PC2')
finish(fig,__file__,'Nearest-neighbor classifiers are fitted in the displayed two-dimensional Iris PCA coordinates. Exact majority votes on a fine grid generate the three decision maps. These are illustrative training-data boundaries, not a held-out accuracy estimate.',data_sources=['data/figure-inputs/iris.npz'],parameters={'R':[1,5,25],'grid_side':200})
