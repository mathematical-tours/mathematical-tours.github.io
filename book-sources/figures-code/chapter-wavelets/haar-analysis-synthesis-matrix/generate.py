from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


N=8;W=np.column_stack([np.concatenate(pywt.dwt(np.eye(N)[i],'haar',mode='periodization')) for i in range(N)])
fig,axs=canvas(1,2,height=3.8)
for ax,A,title in [(axs[0,0],W,r'Analysis: $c=Wf$'),(axs[0,1],W.T,r'Synthesis: $f=W^{\mathsf{T}}c$')]:
    ax.imshow(A,cmap=DIVERGING,vmin=-1/np.sqrt(2),vmax=1/np.sqrt(2));ax.set(title=title,xticks=range(N),yticks=range(N),xlabel='Input coordinate',ylabel='Output coordinate')
    for i,j in zip(*np.nonzero(A)):ax.text(j,i,'+' if A[i,j]>0 else '−',ha='center',va='center',color='white',fontsize=12)
err=np.linalg.norm(W.T@W-np.eye(N));assert err<1e-12
finish(fig,__file__,'Displays the actual N=8 one-step Haar analysis matrix and its transpose. Nonzero entries are signed 1/sqrt(2), and W^T W=I is checked. The layout makes explicit that an orthogonal inverse is the transpose and reverses the order of composed stages.',parameters={'N':N},checks={'orthogonality_residual':float(err)},arrays={'analysis_matrix':W})
