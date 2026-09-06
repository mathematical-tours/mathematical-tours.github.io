from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


from scipy.stats import poisson
k=np.arange(46);fig,axs=canvas(height=3.7,width=7)
for lam in [1,4,10,20]:axs[0,0].plot(k,poisson.pmf(k,lam),'-o',ms=3,label=rf'$\lambda={lam}$')
axs[0,0].set(title='Poisson probability masses',xlabel='Photon count k',ylabel=r'$\mathbb{P}(Y=k)$');axs[0,0].legend();axs[0,0].grid(True)
finish(fig,__file__,'Displays discrete Poisson probability masses, not continuous density curves, for means 1,4,10,20. The parameter is both the expectation and variance. The plotted support is truncated at 45 counts; omitted probability is recorded.',parameters={'means':[1,4,10,20],'maximum_displayed_count':45},checks={'omitted_tail_mass':[float(poisson.sf(45,m)) for m in [1,4,10,20]]})
