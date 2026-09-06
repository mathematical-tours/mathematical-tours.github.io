from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


f,h,H,y=spikes_case();A=lambda x:convolution(x,H);At=lambda x:convolution(x,H.conj());lam=.003;ref=ista(y,A,At,lam,steps=15000,accelerated=True);x,E,it=ista(y,A,At,lam,steps=3000,history=True);Eref=.5*np.linalg.norm(A(ref)-y)**2+lam*np.sum(abs(ref));gap=np.maximum(E-Eref,1e-16);err=np.linalg.norm(it-ref,axis=1)/np.linalg.norm(f);fig,a=canvas(1,2,height=3.8)
curve(a[0,0],np.arange(1,len(E)+1),np.log10(gap/Eref),'Relative objective gap','ISTA iteration',r'$\log_{10}((E_k-E_{ref})/E_{ref})$');curve(a[0,1],np.arange(1,len(E)+1),np.log10(np.maximum(err,1e-16)),'Distance to numerical reference','ISTA iteration',r'$\log_{10}(\|x_k-x_{ref}\|/\|f_0\|)$');assert np.max(np.diff(E))<1e-10
res=np.linalg.norm(ref-soft(ref-At(A(ref)-y),lam));finish(fig,__file__,'Actual ISTA objective values and iterates are compared with a much longer accelerated solve of the same objective. The reference is labeled numerical rather than exact. The unit step is valid because max|H|²=1; monotonicity and the reference proximal fixed-point residual are recorded.',parameters={'lambda':lam,'iterations':3000,'reference_iterations':15000},checks={'maximum_objective_increase':float(np.max(np.diff(E))),'reference_fixed_point_residual':float(res)},arrays={'objective':E,'relative_iterate_error':err})
