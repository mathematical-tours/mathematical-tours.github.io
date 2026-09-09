from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from related_transform_data import notch
x=np.linspace(-1.45,1.45,121);X,Y=np.meshgrid(x,x);Z=abs(notch(X+1j*Y));limit=2
fig=plt.figure(figsize=(6.5,4.9));ax=fig.add_subplot(projection='3d')
surf=ax.plot_surface(X,Y,np.minimum(Z,limit),cmap='viridis',vmin=0,vmax=limit,rcount=90,ccount=90,linewidth=.16,edgecolor='#344b5550',antialiased=True)
t=np.linspace(0,2*np.pi,400);ax.plot(np.cos(t),np.sin(t),np.zeros_like(t),color=TEAL,lw=1.3)
ax.set(xlabel=r'$\operatorname{Re}z$',ylabel=r'$\operatorname{Im}z$',zlabel='$|H(z)|$',zlim=(0,limit),xticks=[-1,0,1],yticks=[-1,0,1]);ax.view_init(26,-124)
ax.set_box_aspect((1,1,.72));fig.colorbar(surf,ax=ax,shrink=.56,pad=.04,label='magnitude (clipped)')
save(fig,'chapter-related-transforms','notch-transfer-surface')
