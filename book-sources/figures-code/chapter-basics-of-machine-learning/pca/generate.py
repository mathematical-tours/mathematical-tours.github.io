from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


X,y=iris();m,s,V,Z=pca(X);fig=plt.figure(figsize=(9,4),layout='constrained');a=fig.add_subplot(121);b=fig.add_subplot(122,projection='3d')
for k,name in enumerate(['setosa','versicolor','virginica']):a.scatter(*Z[y==k,:2].T,s=19,color=COLORS[k],label=name);b.scatter(*Z[y==k,:3].T,s=13,color=COLORS[k])
a.set(title='Two principal components',xlabel='PC1',ylabel='PC2');a.legend();b.set(title='Three principal components',xlabel='PC1',ylabel='PC2',zlabel='PC3')
finish(fig,__file__,'Both views use the same centered Iris data and orthonormal PCA axes, without standardizing the centimeter measurements. Species labels are used only for color, not to fit the projection.',data_sources=['data/figure-inputs/iris.npz'],checks={'two_component_variance_fraction':float(np.sum(s[:2]**2)/np.sum(s*s))})
