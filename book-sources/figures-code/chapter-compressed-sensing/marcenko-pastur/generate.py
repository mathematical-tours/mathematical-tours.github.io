from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


rng=np.random.default_rng(20260906);fig,a=canvas(2,2,height=6);arrays={}
for ax,beta in zip(a.ravel(),[.5,.25,.1,.02]):
 s=160;P=int(s/beta);B=rng.standard_normal((P,s))/np.sqrt(P);ev=np.linalg.eigvalsh(B.T@B);lo=(1-np.sqrt(beta))**2;hi=(1+np.sqrt(beta))**2;x=np.linspace(lo,hi,800);density=np.sqrt(np.maximum((hi-x)*(x-lo),0))/(2*np.pi*beta*x);ax.hist(ev,bins=25,density=True,color=TEAL,alpha=.3,label='empirical');ax.plot(x,density,color=RED,label='MP density');ax.set(title=rf'$\beta=s/P={beta:g}$, s={s}, P={P}',xlabel='Gram eigenvalue',ylabel='Density');ax.legend();arrays[str(beta)]=ev
finish(fig,__file__,'Each histogram comes from eigenvalues of BᵀB with independent N(0,1/P) entries and the stated aspect ratio s/P. The overlay is the exact Marchenko–Pastur density for beta<1, so there is no atom at zero. Density normalization is not confused with singular-value density.',parameters={'seed':20260906,'s':160,'aspect_ratios':[.5,.25,.1,.02]},arrays=arrays)
