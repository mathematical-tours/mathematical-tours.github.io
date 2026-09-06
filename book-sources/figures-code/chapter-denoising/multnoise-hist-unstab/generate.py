from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


from scipy.stats import gamma,norm
K=4;rng=np.random.default_rng(20260906);w=rng.gamma(K,1/K,200000);c=special.digamma(K)-np.log(K);z=np.log(w)-c;variance=float(special.polygamma(1,K));fig,axs=canvas(1,2,height=3.8)
axs[0,0].hist(w,bins=80,density=True,color=TEAL,alpha=.35);x=np.linspace(.001,4,600);axs[0,0].plot(x,gamma.pdf(x,K,scale=1/K),color=TEAL);axs[0,0].set(title='Mean-one multiplicative noise',xlabel=r'$w$',ylabel='Density',xlim=(0,4))
axs[0,1].hist(z,bins=90,density=True,color=TEAL,alpha=.35);x=np.linspace(-3,2,600);v=np.exp(x+c);axs[0,1].plot(x,gamma.pdf(v,K,scale=1/K)*v,label='exact log-Gamma',color=TEAL);axs[0,1].plot(x,norm.pdf(x,scale=np.sqrt(variance)),'--',label='Gaussian approximation',color=RED);axs[0,1].set(title='Centered additive log noise',xlabel=r'$z=\log w-c$',ylabel='Density',xlim=(-3,2));axs[0,1].legend()
finish(fig,__file__,'The histograms use the same Gamma multipliers before and after logarithmic stabilization. Centering uses c=digamma(K)-log(K), and the exact transformed density includes the exponential Jacobian. A Gaussian with variance trigamma(K) is shown only as an approximation; log-Gamma noise is not exactly Gaussian.',parameters={'seed':20260906,'looks':K,'samples':len(w)},checks={'theoretical_log_mean_before_centering':float(c),'theoretical_log_variance':variance,'empirical_centered_mean':float(z.mean()),'empirical_variance':float(z.var())})
