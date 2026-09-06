from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


rng=np.random.default_rng(20260906);base=rng.standard_normal((160,2));labels=np.repeat([0,1],80);fig,a=canvas(1,3,height=3.8)
for ax,separation in zip(a.ravel(),[.5,1.5,3]):
 X=base.copy();X[:,0]+=(2*labels-1)*separation/2;model=LogisticRegression(C=10,max_iter=1000).fit(X,labels);xx,yy,grid=grid2(X,140,.5);prob=model.predict_proba(grid)[:,1].reshape(xx.shape);o=ax.pcolormesh(xx,yy,prob,cmap=DIVERGING,vmin=0,vmax=1,shading='nearest');ax.contour(xx,yy,prob,levels=[.5],colors=INK,linewidths=1)
 for k in [0,1]:ax.scatter(*X[labels==k].T,color=COLORS[k],s=8,edgecolor='white',linewidth=.2)
 ax.set(title=f'Mean separation {separation}',xlabel='x₁',ylabel='x₂')
fig.colorbar(o,ax=a.ravel().tolist(),label='Fitted P(Y=1 | x)',shrink=.7);finish(fig,__file__,'Three logistic models are fitted to paired Gaussian clouds with a shared random realization and increasing class-mean separation. Probability maps use the same 0–1 color scale and finite quadratic regularization, so complete separation does not send an unregularized coefficient to infinity.',parameters={'seed':20260906,'samples':160,'separations':[.5,1.5,3],'sklearn_C':10})
