"""Finite-sum logistic experiments with a checked numerical reference."""
from learning_models import *
from scipy.optimize import minimize

def optimization_runs(epochs=60):
    X,y,lam,loss,grad,grad_i=logistic_case();n,p=X.shape
    fit=minimize(loss,np.zeros(p),jac=grad,method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-12,'maxiter':3000,'maxls':50});ref=fit.x;fref=loss(ref);L=np.linalg.norm(X,2)**2/(4*n)+lam;Li=np.max(np.sum(X*X,axis=1))/4+lam
    rng=np.random.default_rng(20260906);indices=rng.integers(n,size=epochs*n);sgd=np.zeros(p);aux=np.zeros(p);average=np.zeros(p);sag=np.zeros(p);table=np.zeros((n,p));g=np.zeros(p);batch=np.zeros(p);runs={k:[loss(np.zeros(p))-fref] for k in ['SGD','SGA','SAG','Batch']}
    for k,i in enumerate(indices):
        tau=1/Li/(1+k/(3*n));sgd-=tau*grad_i(sgd,i)
        aux-=1/Li/(1+np.sqrt(k/n))*grad_i(aux,i);average+=(aux-average)/(k+1)
        h=grad_i(sag,i);g+=(h-table[i])/n;table[i]=h;sag-=g/(16*Li)
        if (k+1)%n==0:
            batch-=grad(batch)/L
            for name,b in [('SGD',sgd),('SGA',average),('SAG',sag),('Batch',batch)]:runs[name].append(max(loss(b)-fref,1e-16))
    return np.arange(epochs+1),runs,{'reference_gradient_norm':float(np.linalg.norm(grad(ref))),'reference_objective':float(fref),'batch_L':float(L),'individual_L':float(Li),'samples':n,'features':p,'quadratic_penalty':lam}
