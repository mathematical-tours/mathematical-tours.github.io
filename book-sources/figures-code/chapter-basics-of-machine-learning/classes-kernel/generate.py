from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


from sklearn.datasets import make_moons
from scipy.linalg import eigh
X,y=make_moons(n_samples=200,noise=.16,random_state=20260906);xx,yy,grid=grid2(X,180,.4);fig,a=canvas(1,3,height=3.8)
for ax,width in zip(a.ravel(),[.15,.5,1.5]):
 K=np.exp(-np.sum((X[:,None]-X[None,:])**2,axis=2)/(2*width**2));ev,U=eigh(K);keep=ev>1e-10;features=U[:,keep]*np.sqrt(ev[keep]);model=LogisticRegression(C=3,fit_intercept=True,max_iter=1500).fit(features,y);Ktest=np.exp(-np.sum((grid[:,None]-X[None,:])**2,axis=2)/(2*width**2));features_test=Ktest@(U[:,keep]/np.sqrt(ev[keep]));score=model.decision_function(features_test).reshape(xx.shape);ax.contourf(xx,yy,score,levels=[-1e6,0,1e6],colors=[TEAL,RED],alpha=.16);ax.contour(xx,yy,score,levels=[0],colors=INK,linewidths=1)
 for k in [0,1]:ax.scatter(*X[y==k].T,s=9,color=COLORS[k])
 ax.set(title=f'Gaussian kernel σ = {width}',xlabel='x₁',ylabel='x₂')
finish(fig,__file__,'A Gaussian-kernel logistic classifier is fitted independently for each bandwidth on one fixed two-moons sample. Decision regions are the computed sign of its score, not hand-drawn curves or uncalibrated probabilities. The sample is synthetic and labeled as such.',parameters={'seed':20260906,'samples':200,'noise':.16,'kernel_widths':[.15,.5,1.5],'C':3})
