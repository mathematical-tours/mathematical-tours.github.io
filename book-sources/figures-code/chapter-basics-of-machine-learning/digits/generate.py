from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


X,y=digits();m,s,V,Z=pca(X);fig=plt.figure(figsize=(10,3.7),layout='constrained');ax=fig.add_subplot(131);montage=np.vstack([np.hstack([np.pad(X[np.flatnonzero(y==k)[0]].reshape(8,8),1) for k in range(j,j+5)]) for j in [0,5]]);image(ax,montage,'Acquired 8 × 8 digit images');b=fig.add_subplot(132);c=fig.add_subplot(133,projection='3d');colors=plt.get_cmap('tab10')
for k in range(10):sel=y==k;b.scatter(*Z[sel,:2].T,s=3,color=colors(k),alpha=.5,label=str(k));c.scatter(*Z[sel,:3].T,s=2,color=colors(k),alpha=.5)
b.set(title='Two PCA coordinates',xlabel='PC1',ylabel='PC2');c.set(title='Three PCA coordinates',xlabel='PC1',ylabel='PC2',zlabel='PC3');b.legend(ncol=5,loc='upper center',fontsize=7,markerscale=2)
finish(fig,__file__,'Uses the real 1797-image optical handwritten-digits dataset bundled with scikit-learn. Pixel intensities are divided by 16, and both projections use one centered PCA fit. Labels color the plots but do not determine the principal components.',data_sources=['data/figure-inputs/digits.npz'],parameters={'samples':len(X),'features':64,'pixel_scale':16})
