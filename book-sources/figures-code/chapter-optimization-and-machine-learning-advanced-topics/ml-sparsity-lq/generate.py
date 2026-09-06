from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


fig,a=canvas(1,6,height=2.8,width=13)
for ax,q in zip(a.ravel(),[.5,1,1.5,2,4,np.inf]):lq_plot(ax,q)
finish(fig,__file__,'Exact radial parameterizations of the two-dimensional ell-q unit sets preserve equal axis scales. The q=.5 set is nonconvex; q=1 is a diamond, q=1.5 and q=4 have intermediate strictly convex boundaries, q=2 is a disk, and q=infinity is a square.',parameters={'q':[.5,1,1.5,2,4,'infinity']})
