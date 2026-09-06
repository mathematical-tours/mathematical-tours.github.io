"""Datasets and fitted models shared by related learning figures."""
from common import *
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

def iris():
    d=np.load(ROOT/'data/figure-inputs/iris.npz');return d['X'],d['y']
def digits():
    d=np.load(ROOT/'data/figure-inputs/digits.npz');return d['X']/16,d['y']
def pca(X):
    mean=X.mean(axis=0);U,s,V=np.linalg.svd(X-mean,full_matrices=False)
    for row in V:
        if row[np.argmax(abs(row))]<0:row*=-1
    return mean,s,V,(X-mean)@V.T

def grid2(Z,size=160,pad=.3):
    lo=Z.min(axis=0)-pad;hi=Z.max(axis=0)+pad;x,y=np.meshgrid(np.linspace(lo[0],hi[0],size),np.linspace(lo[1],hi[1],size));return x,y,np.c_[x.ravel(),y.ravel()]

def logistic_case():
    rng=np.random.default_rng(20260906);n=240;p=12;X=rng.standard_normal((n,p));beta=rng.standard_normal(p);y=np.where(rng.random(n)<special.expit(X@beta),1.,-1.);lam=.02
    def loss(b):return np.logaddexp(0,-y*(X@b)).mean()+lam/2*np.sum(b*b)
    def grad(b):return X.T@(-y*special.expit(-y*(X@b)))/n+lam*b
    def grad_i(b,i):return -y[i]*special.expit(-y[i]*(X[i]@b))*X[i]+lam*b
    return X,y,lam,loss,grad,grad_i
