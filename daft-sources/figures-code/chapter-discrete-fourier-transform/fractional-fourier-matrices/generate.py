from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from matplotlib.colors import PowerNorm
from spectral_methods import fractional_dft
n=64;eye=np.eye(n)
fig,axes=plt.subplots(2,4,figsize=(9.6,4.6))
for ax,a in zip(axes.flat,[0,.3,.5,.7,.9,1,1.5,2]):
    matrix=np.column_stack([fractional_dft(e,a) for e in eye])
    assert np.allclose(matrix.conj().T@matrix,eye,atol=1e-12)
    if a==1:assert np.allclose(abs(matrix),1/np.sqrt(n))
    im=ax.imshow(abs(matrix),cmap='magma',norm=PowerNorm(.45,vmin=0,vmax=1),interpolation='nearest')
    axes_box(ax);ax.set(title=rf'$\alpha={a:g}$',xticks=[0,n-1],yticks=[0,n-1])
fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.75,pad=.02,label=r'$|(U^\alpha)_{jk}|$',ticks=[0,.125,.5,1])
save(fig,'chapter-discrete-fourier-transform','fractional-fourier-matrices')
