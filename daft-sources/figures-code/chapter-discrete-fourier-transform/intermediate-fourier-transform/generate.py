from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import fractional_dft
n=17;k=np.arange(-(n//2),n//2+1);f=np.exp(-.5*((k-3)/1.15)**2);x=np.fft.ifftshift(f)
fig,axes=plt.subplots(1,5,figsize=(11,2.6),sharey=True)
for ax,a in zip(axes,np.linspace(0,1,5)):
    y=np.fft.fftshift(fractional_dft(x,a));assert np.isclose(np.linalg.norm(y),np.linalg.norm(x))
    ax.plot(k,y.real,'o-',color=TEAL,ms=2,label='real');ax.plot(k,y.imag,'s--',color=RED,ms=2,label='imaginary')
    axes_box(ax,'index');ax.set(title=rf'$\alpha={a:g}$',xticks=[-8,0,8],ylim=(-.6,1.08))
axes[0].legend(loc='upper left',fontsize=8,frameon=False)
save(fig,'chapter-discrete-fourier-transform','intermediate-fourier-transform')
