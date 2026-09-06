from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


z=np.array([.32,.46,.75]);sgn=np.array([1,1,-1]);sigma=.065;x=np.linspace(0,1,1800)
def kernel(x,z):
 d=np.asarray(x)[:,None]-np.asarray(z)[None,:];C=np.exp(-d*d/(4*sigma*sigma));D2=d/(2*sigma*sigma)*C;D1=-D2;D12=(1/(2*sigma*sigma)-d*d/(4*sigma**4))*C;return C,D1,D2,D12
C,D1,D2,D12=kernel(z,z);b=np.linalg.solve(C,sgn);G=np.block([[C,D2],[D1,D12]]);bc=np.linalg.solve(G,np.r_[sgn,np.zeros(3)]);K,_,K2,_=kernel(x,z);F=K@b;V=np.c_[K,K2]@bc;fig,a=canvas(1,2,height=3.6)
for ax,val,title in zip(a.ravel(),[F,V],[r'Fuchs $\widetilde\eta_F$',r'Vanishing derivative $\widetilde\eta_V$']):curve(ax,x,val,title,'Position','Certificate value');ax.scatter(z,sgn,color=RED,zorder=4);ax.axhline(1,color=MUTED,ls='--');ax.axhline(-1,color=MUTED,ls='--');ax.set_ylim(-1.35,1.5)
res=G@bc-np.r_[sgn,np.zeros(3)];assert abs(res).max()<1e-9
finish(fig,__file__,'Gaussian convolution gives C(x,z)=exp(-(x-z)²/(4 sigma²)), up to an irrelevant common normalization. Both precertificates are obtained from the kernel systems in the text, including the mixed derivative block with its correct sign. Derivative constraints enforce stationarity at support points; feasibility is measured separately and is not assumed.',parameters={'locations':z.tolist(),'signs':sgn.tolist(),'sigma':sigma},checks={'value_and_derivative_constraint_error':float(abs(res).max()),'sampled_F_supnorm':float(abs(F).max()),'sampled_V_supnorm':float(abs(V).max())})
