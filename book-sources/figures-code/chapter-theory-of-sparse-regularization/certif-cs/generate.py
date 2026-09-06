from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


rng=np.random.default_rng(20260906);G=rng.standard_normal((48,64));I=[7,19,38,55];signs=np.array([1,-1,1,1]);fig,a=canvas(1,3,height=3.4);records={}
for ax,P in zip(a.ravel(),[8,20,48]):
 A=G[:P]/np.sqrt(P);B=A[:,I];eta=A.T@B@np.linalg.solve(B.T@B,signs);ax.stem(np.arange(64),eta,linefmt=TEAL,markerfmt=' ',basefmt=' ');ax.scatter(I,eta[I],color=RED,s=24);ax.axhline(1,color=MUTED,ls='--');ax.axhline(-1,color=MUTED,ls='--');ax.set(title=f'P = {P}, N = 64',xlabel='Coefficient index',ylim=(-2.6,2.6));records[str(P)]={'off_support_max':float(np.max(abs(np.delete(eta,I)))),'interpolation_error':float(np.max(abs(eta[I]-signs)))}
finish(fig,__file__,'The Fuchs precertificate is recomputed by the exact Gram solve Aᵀ A_I(A_Iᵀ A_I)^(-1)sign(x_I). The same four signed support entries and a nested Gaussian array are used as P increases. Dashed ±1 lines show dual feasibility; no off-support clipping is applied.',parameters={'seed':20260906,'N':64,'P':[8,20,48],'support':I,'signs':signs.tolist()},checks=records)
