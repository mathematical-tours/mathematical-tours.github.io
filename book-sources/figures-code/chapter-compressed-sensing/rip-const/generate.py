from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


rng=np.random.default_rng(20260906);N=64;ks=np.arange(1,21);fig,a=canvas(height=4,width=7);arrays={'k':ks}
for P in [24,40,56]:
 A=rng.standard_normal((P,N))/np.sqrt(P);bounds=[]
 for k in ks:
  vals=[]
  for _ in range(160):I=rng.choice(N,k,replace=False);ev=np.linalg.eigvalsh(A[:,I].T@A[:,I]);vals.append(max(abs(ev-1)))
  bounds.append(max(vals))
 # delta_k is monotone, so previous sampled bounds also bound all larger k.
 bounds=np.maximum.accumulate(bounds);curve(a[0,0],ks,bounds,'Sampled lower bounds on restricted isometry','Support size k',r'$\widehat\delta_k$',label=f'P={P}');arrays[str(P)]=bounds
a[0,0].axhline(1,color=MUTED,ls='--',lw=.8);a[0,0].legend();finish(fig,__file__,'For each support size, the largest observed ||A_IᵀA_I-I||2 is computed over 160 sampled subsets. Cumulative maxima remain rigorous lower bounds by monotonicity of the true RIP constant. These sampled bounds cannot certify a small RIP constant and are not presented as exact values.',parameters={'seed':20260906,'N':64,'P':[24,40,56],'sampled_supports_per_k':160,'entry_variance':'1/P'},arrays=arrays)
