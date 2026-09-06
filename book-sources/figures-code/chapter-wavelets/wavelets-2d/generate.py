from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


phi,psi,z=pywt.Wavelet('db3').wavefun(level=8)
u=np.linspace(-3,3,85);p=np.interp(u+2,z,phi,left=0,right=0);q=np.interp(u+2,z,psi,left=0,right=0)
X,Y=np.meshgrid(u,u);values=[np.outer(p,q),np.outer(q,p),np.outer(q,q)]
fig=plt.figure(figsize=(12,3.5),layout='constrained');titles=[r'$\psi(x_1)\phi(x_2)$',r'$\phi(x_1)\psi(x_2)$',r'$\psi(x_1)\psi(x_2)$']
lim=max(float(abs(v).max()) for v in values)
for k,(v,title) in enumerate(zip(values,titles),1):
    ax=fig.add_subplot(1,4,k,projection='3d')
    ax.plot_surface(X,Y,v,rstride=3,cstride=3,cmap=DIVERGING,vmin=-lim,vmax=lim,edgecolor=MUTED,linewidth=.18,antialiased=True)
    ax.set(title=title,xlabel=r'$x_1$',ylabel=r'$x_2$',zlim=(-lim,lim),xticks=[-3,0,3],yticks=[-3,0,3],zticks=[-1,0,1])
    ax.view_init(elev=27,azim=-56);ax.set_box_aspect((1,1,.65));ax.tick_params(labelsize=7,pad=0)
ax=fig.add_subplot(1,4,4);diagram_axis(ax)
ax.add_patch(Rectangle((.2,.2),.6,.6,facecolor='#E9F2F3',edgecolor=TEAL))
ax.plot(.5,.5,'o',color=TEAL);ax.text(.5,.43,r'$2^j n$',ha='center')
arrow(ax,(.2,.08),(.8,.08),r'$K2^j$',style='<->');ax.text(.5,.95,'Support enclosure',ha='center',fontsize=11)
finish(fig,__file__,'Restores the three surface plots and the separate support enclosure of the original. These are actual tensor products of the Daubechies-3 scaling and wavelet functions, with shared height and spatial scales. An integer shift puts the generator support in [-2,3], hence in the valid centered K=6 enclosure; the enclosure is not claimed to be minimal.',parameters={'wavelet':'db3','integer_generator_shift':2,'support_enclosure_K':6,'surface_grid':85})
