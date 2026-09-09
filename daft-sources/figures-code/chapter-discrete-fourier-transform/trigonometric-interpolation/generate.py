from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import trigonometric_interpolate
x=np.random.default_rng(42).uniform(-1,1,11);eta=100;y=trigonometric_interpolate(x,eta)
assert np.allclose(y[::eta],x)
fig,ax=plt.subplots(figsize=(8,2.9))
ax.plot(np.arange(12),np.r_[x,x[0]],color=MUTED,lw=.8,alpha=.6,label='sample polygon')
ax.plot(np.arange(11*eta+1)/eta,np.r_[y,y[0]],color=TEAL,label='trigonometric interpolant')
ax.scatter(np.arange(11),x,s=22,color=RED,zorder=3,label='samples')
axes_box(ax,'sample coordinate');ax.set_xlim(0,11);ax.legend(frameon=False,ncol=3,fontsize=9,loc='lower left')
save(fig,'chapter-discrete-fourier-transform','trigonometric-interpolation')
