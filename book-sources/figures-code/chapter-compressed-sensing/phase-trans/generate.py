from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


from scipy.optimize import linprog
rng=np.random.default_rng(20260906);N=48

def trial(P,s):
 A=rng.standard_normal((P,N));A/=np.linalg.norm(A,axis=0);I=rng.choice(N,s,replace=False);J=np.setdiff1d(np.arange(N),I);sgn=rng.choice([-1.,1.],s);B=A[:,I];G=B.T@B;cross=B.T@A[:,J];U=np.linalg.solve(G,cross);eta=sgn@U;den=1-np.max(np.sum(abs(G-np.eye(s)),axis=0));weak=den>0 and np.max(np.sum(abs(cross),axis=0))<den;erc=np.max(np.sum(abs(U),axis=0))<1;fuchs=np.max(abs(eta))<=1+1e-9;x=np.zeros(N);x[I]=sgn;sol=linprog(np.ones(2*N),A_eq=np.c_[A,-A],b_eq=A@x,bounds=(0,None),method='highs');success=sol.success and np.linalg.norm(sol.x[:N]-sol.x[N:]-x)<1e-5
 return np.array([weak,erc,fuchs,success],float)
Ps=np.arange(8,41,4);ks=np.arange(1,21);prob=np.zeros((len(Ps),len(ks)))
for p,P in enumerate(Ps):
 for j,k in enumerate(ks):prob[p,j]=np.mean([trial(P,k)[3] for _ in range(12)]) if k<P else 0
ks2=np.arange(1,17);criteria=np.array([np.mean([trial(24,k) for _ in range(60)],axis=0) for k in ks2]);fig,a=canvas(1,2,height=4);o=a[0,0].imshow(prob,origin='lower',extent=[.5,20.5,6,42],aspect='auto',vmin=0,vmax=1,cmap='viridis');a[0,0].set(title='Empirical exact recovery (N=48)',xlabel='Support size s',ylabel='Measurements P');fig.colorbar(o,ax=a[0,0],label='Success probability')
for j,(name,color) in enumerate([('Weak ERC',BLUE),('ERC',INK),('Fuchs',TEAL),('Exact recovery',RED)]):curve(a[0,1],ks2,criteria[:,j],'Recovery criteria at P=24','Support size s','Empirical probability',label=name,color=color)
a[0,1].legend();a[0,1].set_ylim(-.03,1.03)
finish(fig,__file__,'All probabilities are Monte Carlo counts from column-normalized Gaussian matrices and independent random supports/signs. Exact recovery is solved by linear programming with relative-scale tolerance 1e-5. ERC tests max_j||A_I^+ a_j||1<1. Weak ERC uses the Neumann-series bound max_j||A_Iᵀa_j||1 / (1-||A_IᵀA_I-I||1)<1, requiring a positive denominator. Fuchs tests the sign-specific off-support certificate. Finite trial counts are stated; these are not asymptotic threshold curves.',parameters={'seed':20260906,'N':N,'heatmap_trials_per_cell':12,'criterion_trials_per_sparsity':60,'criterion_P':24},arrays={'P':Ps,'s':ks,'exact_probability':prob,'criteria_s':ks2,'criteria_probability':criteria})
