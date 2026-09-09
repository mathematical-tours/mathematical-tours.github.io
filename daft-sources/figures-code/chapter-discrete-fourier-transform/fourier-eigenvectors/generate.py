from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import fourier_eigenbasis
n=128;V=np.fft.fftshift(fourier_eigenbasis(n),axes=0);k=np.arange(-n//2,n//2)
fig,ax=plt.subplots(figsize=(8,3))
for j,color in enumerate([TEAL,RED,BLUE]):ax.plot(k,abs(V[:,j]),color=color,label=rf'$|v_{j}|$')
axes_box(ax,'centered sample index','magnitude');ax.set(xlim=(-28,28),ylim=(0,.46));ax.legend(frameon=False,ncol=3)
save(fig,'chapter-discrete-fourier-transform','fourier-eigenvectors')
