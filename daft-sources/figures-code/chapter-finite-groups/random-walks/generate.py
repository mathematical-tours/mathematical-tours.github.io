"""Two probability kernels on Z/15Z, with the same four observation times."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'tools'))
from figure_style import *

n=15
v1=np.zeros(n); v1[[1,-1]]=.5
v2=np.zeros(n); v2[[0,1,2,-2,-1]]=[.3,.3,.2,.1,.1]
fig,axs=plt.subplots(4,2,figsize=(7.2,5.0),sharex=True,sharey='row')
for col,v in enumerate([v1,v2]):
    transition=np.array([[v[(j-i)%n] for i in range(n)] for j in range(n)])
    for row,k in enumerate([0,1,5,20]):
        p=np.fft.ifft(np.fft.fft(v)**k).real
        assert np.allclose(p,np.linalg.matrix_power(transition,k)[:,0])
        assert np.isclose(p.sum(),1) and p.min()>-1e-12
        ax=axs[row,col]
        discrete(ax,np.arange(n),p,TEAL if col==0 else BLUE)
        ax.axhline(1/n,color=RED,lw=.8,ls='--')
        ax.text(.5,.82,rf'$k={k}$',ha='center',transform=ax.transAxes)
        ax.set(xlim=(-.4,14.4),xticks=[0,5,10,14])
        axes_box(ax,'$x$' if row==3 else None,'$p^{(k)}(x)$' if col==0 else None)
        if row==0:
            ax.set_title('Nearest-neighbor walk' if col==0 else 'Asymmetric walk')
axs[3,1].text(.98,.55,'Uniform law',color=RED,fontsize=8,ha='right',transform=axs[3,1].transAxes)
save(fig,'chapter-finite-groups','random-walks')
