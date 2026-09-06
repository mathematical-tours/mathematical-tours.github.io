from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


fig,axs=canvas(1,4,height=3.1);ax=axs[0,0];m=np.arange(-4,5);xx,yy=np.meshgrid(m,m);ax.scatter(xx,yy,s=7,color=MUTED)
frequencies=[(1,0),(0,2),(2,3)]
for col,(kx,ky) in zip(COLORS,frequencies):ax.scatter([kx],[ky],color=col,s=45);ax.text(kx+.15,ky+.18,f'({kx},{ky})',fontsize=9,color=col)
ax.set(title='Frequency lattice',xlabel=r'$k_1$',ylabel=r'$k_2$',aspect='equal');ax.axhline(0,lw=.5,color=MUTED);ax.axvline(0,lw=.5,color=MUTED)
y,x=np.mgrid[:192,:192]/192
for ax,(kx,ky) in zip(axs.ravel()[1:],frequencies):image(ax,np.cos(2*np.pi*(kx*x+ky*y)),rf'$\cos(2\pi({kx}x_1+{ky}x_2))$',signed=True)
finish(fig,__file__,'The complex Fourier basis is indexed by integer lattice points. The image panels show exact real parts cos(2 pi (k1 x1+k2 x2)); their orientation and oscillation count follow the marked lattice vectors. All panels have the same spatial domain and amplitude range.',parameters={'frequencies':frequencies,'grid':192})
