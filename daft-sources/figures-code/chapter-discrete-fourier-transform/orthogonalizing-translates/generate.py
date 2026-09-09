from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import orthogonalize_translates
n=20;k=np.arange(-n//2,n//2)
fig,axes=plt.subplots(2,2,figsize=(8.5,4.6),sharex=True)
for row,a in enumerate([.5,.9]):
    f=np.zeros(n);f[0]=1;f[1]=a;g=orthogonalize_translates(f)
    assert np.allclose(abs(np.fft.fft(g)),np.sqrt(n))
    for ax,y,title in zip(axes[row],[f,g],[rf'Input: $a={a}$',rf'Orthogonalized: $a={a}$']):
        discrete(ax,k,np.fft.fftshift(y));axes_box(ax,'centered index');ax.set_title(title);ax.axhline(0,color=MUTED,lw=.5)
    axes[row,0].set_ylim(-.1,1.15);axes[row,1].set_ylim(-1.3,4.5)
save(fig,'chapter-discrete-fourier-transform','orthogonalizing-translates')
