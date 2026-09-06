from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


rng=np.random.default_rng(20260906);N=2**np.arange(4,17);trials=384;samples=np.abs(rng.standard_normal((trials,int(N[-1]))));maxima=np.maximum.accumulate(samples,axis=1)[:,N-1];mean=maxima.mean(axis=0);std=maxima.std(axis=0,ddof=1);fig,axs=canvas(2,1,height=5.4,width=7)
curve(axs[0,0],N,mean,'Largest absolute Gaussian coefficient','N','Empirical mean',label='Monte Carlo');axs[0,0].plot(N,np.sqrt(2*np.log(N)),'--',label=r'$\sqrt{2\log N}$',color=RED);axs[0,0].legend();curve(axs[1,0],N,std,'Fluctuations shrink as N increases','N','Empirical standard deviation');
for ax in axs.ravel():ax.set_xscale('log',base=2)
finish(fig,__file__,'For each N, computes max_i |Z_i| over 384 independent standard-Gaussian vectors and reports its empirical mean and standard deviation. Prefixes share draws across N only to smooth the comparison; trials remain independent. The dashed leading-order scale is sqrt(2 log N). This is the absolute maximum, unlike the signed maximum in the compressed-sensing chapter.',parameters={'seed':20260906,'trials':trials,'dimensions':N.tolist()},checks={'all_standard_deviations_nonnegative':bool(np.all(std>=0))},arrays={'dimensions':N,'maxima':maxima,'mean':mean,'std':std})
