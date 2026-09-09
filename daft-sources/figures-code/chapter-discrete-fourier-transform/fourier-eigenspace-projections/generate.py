from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import fourier_components
n=17;k=np.arange(-(n//2),n//2+1);f=np.exp(-.5*((k-3)/1.15)**2);x=np.fft.ifftshift(f)
components=np.fft.fftshift(fourier_components(x),axes=1)
assert np.allclose(components.sum(axis=0),f)
assert np.max(abs(components.imag))<1e-14
fig,axes=plt.subplots(2,2,figsize=(8,4.5),sharex=True,sharey=True)
for ax,y,label in zip(axes.flat,components,[r'$P_1 f$',r'$P_{-1} f$',r'$P_{i} f$',r'$P_{-i} f$']):
    discrete(ax,k,y.real);ax.axhline(0,color=MUTED,lw=.5);axes_box(ax,'centered index');ax.set_title(label)
save(fig,'chapter-discrete-fourier-transform','fourier-eigenspace-projections')
