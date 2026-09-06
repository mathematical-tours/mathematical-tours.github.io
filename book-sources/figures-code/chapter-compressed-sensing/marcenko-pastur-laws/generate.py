from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from sparse_models import *


fig,a=canvas(height=4,width=7)
for beta in [.1,.25,.5,.8]:lo=(1-np.sqrt(beta))**2;hi=(1+np.sqrt(beta))**2;x=np.linspace(lo,hi,2000);d=np.sqrt(np.maximum((hi-x)*(x-lo),0))/(2*np.pi*beta*x);curve(a[0,0],x,d,'Marchenko–Pastur probability densities','Eigenvalue λ',r'$f_\beta(\lambda)$',label=rf'$\beta={beta}$')
a[0,0].legend();finish(fig,__file__,'Exact Marchenko–Pastur densities are evaluated only on their support [(1-sqrt(beta))²,(1+sqrt(beta))²]. All displayed beta values are below one, avoiding an unrepresented point mass at zero.',parameters={'aspect_ratios':[.1,.25,.5,.8]})
