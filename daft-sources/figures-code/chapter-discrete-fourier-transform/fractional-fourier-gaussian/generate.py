from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import fractional_dft
n=65;k=np.arange(-(n//2),n//2+1);x=np.fft.ifftshift(np.exp(-.5*(k/6)**2))
fig,axes=plt.subplots(3,3,figsize=(9,6),sharex=True,sharey=True)
for ax,a in zip(axes.flat,np.linspace(0,2,9)):
    y=np.fft.fftshift(fractional_dft(x,a))
    ax.plot(k,y.real,color=TEAL,label='real');ax.plot(k,y.imag,'--',color=RED,label='imaginary')
    axes_box(ax,'centered index');ax.set(title=rf'$\alpha={a:g}$',xlim=(-24,24),ylim=(-.6,1.95))
    assert np.isclose(np.linalg.norm(y),np.linalg.norm(x))
axes[0,0].legend(frameon=False,fontsize=8)
save(fig,'chapter-discrete-fourier-transform','fractional-fourier-gaussian')
