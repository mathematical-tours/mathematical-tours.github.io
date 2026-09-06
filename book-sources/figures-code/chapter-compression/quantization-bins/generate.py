from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


u=np.linspace(-3.5,3.5,1401);q=np.sign(u)*np.floor(abs(u));v=np.sign(q)*(abs(q)+.5)
fig,axs=canvas(1,2,height=3.5)
for ax,y,title in [(axs[0,0],q,r'Integer symbol $Q_T(u)$'),(axs[0,1],v,r'Reconstruction $D_T(Q_T(u))$')]:
    for k in range(len(u)-1):
        if y[k]==y[k+1]:ax.plot(u[k:k+2],y[k:k+2],color=TEAL,lw=1.7)
    ax.set(title=title,xlabel=r'$u/T$',xticks=range(-3,4));ax.axvspan(-1,1,color=TEAL,alpha=.08);ax.grid(True)
    for k in [-3,-2,-1,1,2,3]:
        exact=np.sign(k)*np.floor(abs(k));val=exact if title.startswith('Integer') else np.sign(exact)*(abs(exact)+.5);inside=val-np.sign(k) if abs(k)>1 or title.startswith('Integer') else 0
        ax.plot(k,val,'o',color=TEAL,ms=4);ax.plot(k,inside,'o',mfc='white',mec=TEAL,ms=4)
finish(fig,__file__,'Unnumbered dead-zone quantization diagram, in units of T. The zero cell is (-T,T). At |u|=T the symbol is nonzero; nonzero cells are reconstructed at their midpoints sign(q)(|q|+1/2)T, while zero is reconstructed exactly as zero. Endpoint dots distinguish the equality convention from strict hard thresholding.')
