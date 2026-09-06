from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from approximation_models import *


f=flower();M=np.unique(np.geomspace(16,32000,60).astype(int));fig,a=canvas(height=4,width=7);records={}
for kind in ['Fourier','DCT','Local DCT','Wavelet']:
    c,inv=transform(f,kind);err=np.linalg.norm(inv(c)-f)/np.linalg.norm(f);assert err<1e-10;p=np.sort(abs(c).ravel())[::-1]**2;e=np.r_[np.cumsum(p[::-1])[::-1],0][M]/np.sum(f*f);curve(a[0,0],M,e,'Equal real-coefficient budgets','Retained coefficients M','Relative squared error',label=kind);records[kind]=e
ax=a[0,0];ax.set(xscale='log',yscale='log');ax.legend();finish(fig,__file__,'Actual best-M errors in four orthonormal bases use one common image and real coefficient counts. The local DCT consists of disjoint 16 by 16 blocks. Full-transform inversion is checked before computing the discarded-energy curves.',data_sources=['data/flower.png'],parameters={'size':256,'local_dct_block':16},arrays={'M':M,**records})
