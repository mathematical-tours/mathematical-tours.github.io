"""A correctly centered even function and its real Fourier transform on Z/17Z."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'tools'))
from figure_style import *

x = np.arange(-8,9)
f = np.exp(-2*(np.pi*x/17)**2)
# ifftshift converts centered group indices to 0,...,16; the group convention
# has a positive exponential, hence n*ifft rather than fft.
fhat = np.fft.fftshift(17*np.fft.ifft(np.fft.ifftshift(f)))
assert np.max(np.abs(fhat.imag)) < 1e-12
assert np.isclose(fhat[8], f.sum())
fig, axs = plt.subplots(1,2,figsize=(7.0,2.55))
for ax,y,title,label in zip(axs,[f,fhat.real],['Function','Fourier coefficients'],['$x$','$k$']):
    discrete(ax,x,y)
    ax.set(title=title,xticks=[-8,-4,0,4,8],xlim=(-8.5,8.5))
    axes_box(ax,label,'$f(x)$' if label=='$x$' else r'$\widehat f(\chi_k)$')
axs[1].annotate(r'$\widehat f(\chi_0)=\sum_x f(x)$', xy=(0,f.sum()),
                xytext=(2,4.8),fontsize=9,arrowprops={'arrowstyle':'-','lw':.7,'color':MUTED})
save(fig,'chapter-finite-groups','fourier-decomposition')
