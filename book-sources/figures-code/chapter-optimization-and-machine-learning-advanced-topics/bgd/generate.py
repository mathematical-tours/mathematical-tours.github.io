from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


from scipy.optimize import minimize
X,y,lam,loss,grad,grad_i=logistic_case();fit=minimize(loss,np.zeros(X.shape[1]),jac=grad,method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-12,'maxiter':3000});fref=loss(fit.x);L=np.linalg.norm(X,2)**2/(4*len(X))+lam;fig,a=canvas(height=4,width=7);arrays={}
for factor in [.25,1,1.5]:
 x=np.zeros(X.shape[1]);values=[]
 for k in range(180):values.append(max(loss(x)-fref,1e-16));x-=factor/L*grad(x)
 curve(a[0,0],np.arange(180),values,'Batch gradient descent: measured objective gap','Full-gradient iterations','Objective gap',label=rf'$\tau={factor}/L$');arrays[str(factor)]=np.array(values)
a[0,0].set_yscale('log');a[0,0].legend();finish(fig,__file__,'Batch gradient descent is evaluated on one regularized logistic finite-sum problem. The quadratic penalty ensures a finite unique minimizer; a high-accuracy independent solve supplies a numerical objective reference. All three steps are below 2/L, and gaps below 1e-16 are displayed at the numerical floor.',parameters={'seed':20260906,'n':len(X),'p':X.shape[1],'lambda':lam,'iterations':180},checks={'reference_gradient_norm':float(np.linalg.norm(grad(fit.x))),'L':float(L)},arrays=arrays)
