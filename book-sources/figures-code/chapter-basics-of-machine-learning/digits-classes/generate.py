from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


from matplotlib.colors import to_rgb
X,y=digits();train,test=train_test_split(np.arange(len(X)),test_size=.25,random_state=20260906,stratify=y);m,s,V,Z=pca(X[train]);model=LogisticRegression(C=3,max_iter=2500).fit(X[train],y[train]);xx,yy,grid=grid2(Z[:,:2],100,.3);lift=m+grid@V[:2];p=model.predict_proba(lift);colors=np.array([plt.get_cmap('tab10')(k)[:3] for k in range(10)]);fig=plt.figure(figsize=(10,5.8),layout='constrained');gs=fig.add_gridspec(4,6,height_ratios=[1,1,1,.09])
for k in range(9):ax=fig.add_subplot(gs[k//3,k%3]);im=ax.imshow(p[:,k].reshape(xx.shape),origin='lower',extent=[xx.min(),xx.max(),yy.min(),yy.max()],vmin=0,vmax=1,cmap='viridis',aspect='auto');ax.set(title=f'P(class {k})',xticks=[],yticks=[])
fig.colorbar(im,cax=fig.add_subplot(gs[3,:3]),orientation='horizontal',ticks=[0,.5,1],label='Class probability (shared scale)')
ax=fig.add_subplot(gs[:3,3:]);ax.imshow((p@colors).reshape(*xx.shape,3),origin='lower',extent=[xx.min(),xx.max(),yy.min(),yy.max()],aspect='auto');ax.set(title='Mixture of all ten class colors',xlabel='PC1',ylabel='PC2')
projected_test=(X[test]-m)@V[:2].T
ax.scatter(projected_test[:,0],projected_test[:,1],c=colors[y[test]],s=7,edgecolors='black',linewidths=.25,alpha=.85)
ax.set_xlim(xx.min(),xx.max());ax.set_ylim(yy.min(),yy.max())
for k in range(10):ax.scatter([],[],color=colors[k],label=str(k))
ax.legend(ncol=5,loc='lower center',fontsize=8);fig.suptitle('Probabilities on the affine two-PC slice through the training mean',fontsize=12)
finish(fig,__file__,'A ten-class multinomial logistic model is trained on 75% of the real digit images. The nine small panels show classes 0–8 on one probability scale; the RGB mixture includes all ten probabilities. Grid points are lifted from the two-PC plane back to 64 pixel coordinates before prediction. Dots are held-out images projected onto these axes and colored by their true class. Their full-space predictions need not coincide with the affine slice at the projected locations. PCA is fit on training data only.',data_sources=['data/figure-inputs/digits.npz'],parameters={'split_seed':20260906,'test_fraction':.25,'C':3,'grid_side':100},checks={'held_out_accuracy':float(model.score(X[test],y[test])),'probability_sum_error':float(np.max(abs(p.sum(axis=1)-1)))})
