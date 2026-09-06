from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


f=flower(128);mask=np.random.default_rng(20260906).random(f.shape)>.5;y=f*mask;fig,a=canvas(1,4,height=3.6);image(a[0,0],y,'Half of the pixels observed');image(a[0,1],inpaint(y,mask,'sobolev',1000),'Sobolev equality constraint');records={}
for ax,stationary in zip(a[0,2:],[False,True]):
 W,S=wavelet_operator(f.shape,stationary);delta=np.zeros_like(f);delta[0,0]=1;weights=np.linalg.norm(W(delta),axis=(1,2),keepdims=True) if stationary else 1;c=ista(y,lambda c:S(c)*mask,lambda x:W(x*mask),.0005*weights,steps=1000,accelerated=True);g=S(c);name='Stationary frame' if stationary else 'Orthogonal wavelets';image(ax,g,f'{name}\n{snr(f,g):.1f} dB');records[name]={'snr_db':snr(f,g),'observed_residual':float(np.linalg.norm((g-y)*mask))}
finish(fig,__file__,'Reconstructs one masked image with Sobolev interpolation and two coefficient-space ell-1 synthesis models. Wavelet solutions use a small positive lambda and hence have a measured nonzero data residual; they are not labeled exact interpolants. The stationary penalty weights coefficients by their atom L2 norms, equivalent to a dictionary of unit-norm atoms. The stationary frame uses normalized analysis and its verified adjoint synthesis.',data_sources=['data/flower.png'],parameters={'size':128,'seed':20260906,'lambda':.0005,'iterations':1000},checks=records)
