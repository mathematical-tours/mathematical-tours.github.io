from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


X,y=iris();m,s,V,Z=pca(X);rng=np.random.default_rng(20260906);centers=X[[0,60,120]].copy();history=[centers.copy()];energies=[]
for _ in range(12):
 labels=np.argmin(np.sum((X[:,None]-centers[None,:])**2,axis=2),axis=1);centers=np.array([X[labels==k].mean(axis=0) if np.any(labels==k) else centers[k] for k in range(3)]);history.append(centers.copy());energies.append(float(np.sum((X-centers[labels])**2)))
labels=np.argmin(np.sum((X[:,None]-centers[None,:])**2,axis=2),axis=1);fig,a=canvas(1,2,height=4)
for k in range(3):
 a[0,0].scatter(*Z[labels==k,:2].T,color=COLORS[k],s=15,alpha=.6);path=(np.array(history)[:,k]-m)@V[:2].T;a[0,0].plot(*path.T,'-o',color=COLORS[k],ms=4);a[0,0].scatter(*path[-1],marker='X',s=100,color=COLORS[k],edgecolor='white')
counts=np.array([[np.sum((labels==k)&(y==j)) for j in range(3)] for k in range(3)]);bottom=np.zeros(3)
for j,name in enumerate(['setosa','versicolor','virginica']):a[0,1].bar(range(3),counts[:,j],bottom=bottom,color=COLORS[j],label=name);bottom+=counts[:,j]
a[0,0].set(title='Centroid trajectories, viewed by PCA',xlabel='PC1',ylabel='PC2');a[0,1].set(title='True species within fitted clusters',xlabel='Cluster',ylabel='Number of observations',xticks=range(3));a[0,1].legend();assert np.max(np.diff(energies))<1e-8
finish(fig,__file__,'Lloyd iterations run in the original four-dimensional feature space; PCA is only the display projection. Right-hand bars count the true species inside the final fitted clusters, preserving the distinction between unsupervised cluster IDs and class labels.',data_sources=['data/figure-inputs/iris.npz'],parameters={'initial_observation_indices':[0,60,120],'iterations':12},checks={'objective_values':energies,'class_counts_per_cluster':counts.tolist()})
