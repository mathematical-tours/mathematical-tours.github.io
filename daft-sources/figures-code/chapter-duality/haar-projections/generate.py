from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
x=(np.arange(2048)+.5)/2048
f=.25+.8*x+.12*np.sin(2*np.pi*x)+.14*np.sin(5*np.pi*x)
fig,axs=plt.subplots(3,2,figsize=(7.1,5.1),sharex=True)
for row,J in enumerate([0,1,2]):
    ax=axs[row,0];size=len(x)//2**J
    p=np.repeat(f.reshape(2**J,size).mean(axis=1),size)
    ax.plot(x,f,color=MUTED,lw=1,label='$f$')
    ax.step(np.r_[0,x,1],np.r_[p[0],p,p[-1]],where='mid',color=TEAL,label=rf'$P_{{F_{J}}}f$')
    ax.set(ylim=(.1,1.25),yticks=[.2,.6,1.0])
    ax.legend(loc='upper left',frameon=False,fontsize=9,ncol=2)
    ax=axs[row,1]
    if J==0:
        ax.plot([0,1],[1,1],color=TEAL,label=r'$\psi_0$')
        ax.set(ylim=(-.15,1.45),yticks=[0,1])
    else:
        j=J-1
        for k,color in zip(range(2**j),[TEAL,RED]):
            left=k/2**j;mid=(k+.5)/2**j;right=(k+1)/2**j
            xx=[0,left,left,mid,mid,right,right,1]
            y=[0,0,2**(j/2),2**(j/2),-2**(j/2),-2**(j/2),0,0]
            ax.plot(xx,y,color=color,label=rf'$\psi_{{{j},{k}}}$')
        ax.set(ylim=(-1.8,2.1),yticks=[-1,0,1])
    ax.legend(loc='upper right',frameon=False,fontsize=9,ncol=2)
    for col in range(2):
        axs[row,col].set(xlim=(0,1),xticks=[0,.5,1]);axes_box(axs[row,col],'$x$' if row==2 else None)
axs[0,0].set_title('Orthogonal projection');axs[0,1].set_title('New basis functions')
save(fig,'chapter-duality','haar-projections')
