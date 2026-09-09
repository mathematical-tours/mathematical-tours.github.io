from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from matplotlib.colors import LogNorm
n=128;i,j=np.indices((n,n));f=((i<n//2)==(j<n//2)).astype(float)
F=abs(np.fft.fftshift(np.fft.fft2(f)))/n**2
assert np.isclose(F[n//2,n//2],.5)
fig,axes=plt.subplots(1,2,figsize=(7.8,3.3))
axes[0].imshow(f,cmap='gray',vmin=0,vmax=1,interpolation='nearest',extent=(-.5,n-.5,n-.5,-.5))
axes[0].set_title('Spatial pattern');axes_box(axes[0],'column','row')
im=axes[1].imshow(np.maximum(F,1e-4),norm=LogNorm(vmin=1e-4,vmax=.5),cmap='magma',extent=(-.5,.5,.5,-.5),interpolation='nearest')
axes[1].set_title('Centered Fourier magnitude');axes_box(axes[1],'horizontal frequency','vertical frequency')
fig.colorbar(im,ax=axes[1],shrink=.8,label=r'$|\widehat f|/N^2$')
save(fig,'chapter-discrete-fourier-transform','two-dimensional-dft')
