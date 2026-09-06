from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


X,y=iris();m,s,V,Z=pca(X);C=(X-m).T@(X-m)/len(X);ev=s*s/len(X);fig,a=canvas(1,2,height=4);o=a[0,0].imshow(C,cmap='viridis');a[0,0].set(title='Empirical covariance (Iris)',xticks=range(4),yticks=range(4),xticklabels=['Sepal L','Sepal W','Petal L','Petal W'],yticklabels=['Sepal L','Sepal W','Petal L','Petal W']);plt.setp(a[0,0].get_xticklabels(),rotation=35,ha='right');fig.colorbar(o,ax=a[0,0],label='cm²');a[0,1].bar(np.arange(1,5),ev,color=TEAL);a[0,1].set(title='Covariance eigenvalues',xlabel='Principal component',ylabel='Variance (cm²)',xticks=range(1,5));assert np.allclose(np.linalg.eigvalsh(C)[::-1],ev)
finish(fig,__file__,'The 150 real Iris observations are centered without feature standardization. Covariance uses the chapter normalization 1/n. The spectrum displays covariance eigenvalues sigma_k²/n for singular values of the unnormalized centered data, avoiding a missing factor n.',data_sources=['data/figure-inputs/iris.npz'],parameters={'n':len(X),'covariance_normalization':'1/n','feature_scaling':'original centimeters'},checks={'eigenvalues':ev.tolist()})
