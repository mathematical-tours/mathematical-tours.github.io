from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


rng=np.random.default_rng(20260906);N=np.unique(np.geomspace(2,32768,32).astype(int));z=rng.standard_normal((512,N[-1]));m=np.maximum.accumulate(z,axis=1)[:,N-1];mean=m.mean(axis=0);err=m.std(axis=0,ddof=1)/np.sqrt(len(m));fig,a=canvas(height=4,width=7);ax=a[0,0];ax.plot(N,mean,label='Empirical mean of signed maximum');ax.fill_between(N,mean-1.96*err,mean+1.96*err,color=TEAL,alpha=.15,label='Approx. 95% Monte Carlo CI');ax.plot(N,np.sqrt(2*np.log(N)),'--',color=RED,label=r'$\sqrt{2\log N}$');ax.set(xscale='log',title='Maximum of N independent Gaussian variables',xlabel='N',ylabel='Maximum');ax.legend();ax.grid(True)
finish(fig,__file__,'Uses the signed maximum max(Z1,...,ZN), as the caption states, rather than the maximum absolute value shown earlier in the denoising chapter. Independent trials determine the empirical mean and Monte Carlo standard error; the dashed expression is an asymptotic comparison, not an exact finite-N mean.',parameters={'seed':20260906,'independent_trials':512},arrays={'N':N,'mean':mean,'mean_standard_error':err})
