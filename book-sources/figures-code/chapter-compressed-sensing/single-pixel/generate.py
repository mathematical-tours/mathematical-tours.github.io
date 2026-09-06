from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import wavelet_operator,ista
n=64;Q=n*n;f=flower(n);W,S=wavelet_operator(f.shape);rng=np.random.default_rng(20260906);order=np.r_[0,rng.permutation(np.arange(1,Q))];scramble=rng.permutation(Q)
def hadamard(v):
 z=np.asarray(v).copy();h=1
 while h<len(z):
  b=z.reshape(-1,2*h);left=b[:,:h].copy();right=b[:,h:].copy();b[:,:h]=left+right;b[:,h:]=left-right;h*=2
 return z/np.sqrt(len(z))
fig,a=canvas(1,4,height=3.5,width=12);image(a[0,0],f,r'Flower: $Q=4096$ pixels');records=[];arrays={'reference':f}
for ax,P in zip(a.ravel()[1:],[512,1024,2048]):
 rows=order[:P]
 def Phi(x):return hadamard(x.ravel()[scramble])[rows]
 def Phit(y):
  z=np.zeros(Q);z[rows]=y;out=np.empty(Q);out[scramble]=hadamard(z);return out.reshape(n,n)
 A=lambda c:Phi(S(c));At=lambda y:W(Phit(y));y=Phi(f)
 c=ista(y,A,At,.0005,L=1,steps=2200,accelerated=True);g=S(c);image(ax,g,f'P = {P} ({100*P/Q:.1f}%)\n{snr(f,g):.1f} dB')
 v=rng.standard_normal(P);x=rng.standard_normal(f.shape);adj=abs(np.dot(Phi(x),v)-np.sum(x*Phit(v)));assert adj<1e-10
 records.append({'P':P,'snr_db':snr(f,g),'relative_data_residual':float(np.linalg.norm(Phi(g)-y)/np.linalg.norm(y)),'adjoint_residual':float(adj)})
 arrays[str(P)]=g
finish(fig,__file__,'Simulates compressed sensing of the supplied flower with three nested measurement counts. The sensing rows are orthonormal Walsh-Hadamard patterns with a fixed random pixel permutation; signed measurements can be formed by complementary exposures. Every reconstruction solves the stated wavelet-synthesis l1 problem from its own measurements using the same lambda and iterations. The operator and adjoint are exact, with norm one. These are simulated acquisitions, not new camera measurements.',parameters={'seed':20260906,'side':n,'Q':Q,'P':[512,1024,2048],'lambda':.0005,'iterations':2200,'wavelet':'db2','levels':3},checks={'runs':records},arrays=arrays)
