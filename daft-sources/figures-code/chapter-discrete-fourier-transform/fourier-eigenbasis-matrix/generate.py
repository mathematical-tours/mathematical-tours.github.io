from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import fourier_eigenbasis
n=128;V=fourier_eigenbasis(n)
assert np.allclose(V.T@V,np.eye(n),atol=1e-12)
fig,ax=plt.subplots(figsize=(6,4.5))
im=ax.imshow(abs(np.fft.fftshift(V,axes=0)),cmap='magma',vmin=0,vmax=.5,extent=(-.5,n-.5,n/2-.5,-n/2-.5),interpolation='nearest')
axes_box(ax,'basis vector index','centered sample index');fig.colorbar(im,ax=ax,label='magnitude',shrink=.82)
save(fig,'chapter-discrete-fourier-transform','fourier-eigenbasis-matrix')
