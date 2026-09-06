from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from learning_models import *


rng=np.random.default_rng(20260906);trials=1500;steps=160;x=np.full(trials,2.);paths=[x.copy()]
for k in range(steps):tau=.25/(1+k/40);x-=tau*(x-rng.standard_normal(trials));paths.append(x.copy())
paths=np.array(paths);edges=np.linspace(-1.1,2.1,110);density=np.array([np.histogram(z,bins=edges,density=True)[0] for z in paths]);fig,a=canvas(1,2,height=4);from matplotlib.colors import LogNorm
o=a[0,0].imshow(np.ma.masked_less_equal(density.T,0),norm=LogNorm(vmin=.02,vmax=8),origin='lower',extent=[0,steps,edges[0],edges[-1]],aspect='auto',cmap='viridis');a[0,0].set(title='Empirical distribution of iterates',xlabel='Iteration k',ylabel='x');fig.colorbar(o,ax=a[0,0],label='Density (log color scale)');a[0,1].plot(paths[:,:18],alpha=.4,lw=.8);a[0,1].plot(paths.mean(axis=1),color=INK,lw=2,label='Mean of 1500 trials');a[0,1].axhline(0,color=MUTED,ls='--');a[0,1].set(title='Independent stochastic trajectories',xlabel='Iteration k',ylabel='x');a[0,1].legend()
finish(fig,__file__,'Simulates SGD for E[(x-Z)²/2] with Z standard Gaussian, whose minimizer is exactly zero. Every step draws fresh independent noise. The left heatmap is the empirical density of the same trials whose first 18 trajectories appear at right; the deterministic starting point is represented by one histogram bin; a logarithmic color scale resolves the later spread without allowing this initial point mass to dominate the range.',parameters={'seed':20260906,'trials':trials,'steps':steps,'initial_x':2,'schedule':'.25/(1+k/40)'},checks={'final_mean':float(x.mean()),'final_variance':float(x.var())})
