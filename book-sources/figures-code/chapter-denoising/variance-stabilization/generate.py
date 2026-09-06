from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


from scipy.stats import poisson
means=np.geomspace(.005,100,180);curves=[];tail=[]
transforms=[('Anscombe',lambda k:2*np.sqrt(k+3/8)),('Freeman–Tukey',lambda k:np.sqrt(k+1)+np.sqrt(k)),('2 sqrt(x)',lambda k:2*np.sqrt(k))]
for name,g in transforms:
    values=[]
    for lam in means:
        k=np.arange(int(np.ceil(lam+15*np.sqrt(lam+1)+45)));p=poisson.pmf(k,lam);z=g(k);m=np.sum(p*z);values.append(np.sum(p*(z-m)**2));tail.append(float(poisson.sf(k[-1],lam)))
    curves.append(values)
fig,axs=canvas(height=3.9,width=7)
for (name,_),values in zip(transforms,curves):curve(axs[0,0],means,values,'Poisson variance stabilization','Mean count λ','Transformed variance',label=name)
axs[0,0].set_xscale('log');axs[0,0].axhline(1,color=MUTED,ls='--',lw=.9);axs[0,0].legend();assert max(tail)<1e-14
finish(fig,__file__,'Computes each variance directly from the Poisson probability mass function, including the transformed expectation rather than replacing it by the transform of the mean. The summation tail is below 1e-14. The horizontal line is the asymptotic unit variance; the plot also resolves the low-count regime where stabilization is imperfect.',parameters={'mean_range':[.005,100],'mean_samples':180},checks={'maximum_omitted_tail':max(tail),'nonnegative_variances':bool(np.min(curves)>=0)},arrays={'mean':means,'variances':np.array(curves)})
