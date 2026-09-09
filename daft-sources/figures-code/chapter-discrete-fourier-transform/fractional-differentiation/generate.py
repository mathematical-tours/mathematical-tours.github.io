from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from spectral_methods import fractional_derivative
n=513;period=24;x=(np.arange(n)-n//2)*period/n;f=np.exp(-x*x/2)
fig,axes=plt.subplots(3,3,figsize=(9,6),sharex=True,sharey=True)
for ax,a in zip(axes.flat,np.linspace(0,2,9)):
    y=fractional_derivative(f,a,period)
    ax.plot(x,y,color=TEAL);ax.axhline(0,color=MUTED,lw=.5);axes_box(ax,'$x$');ax.set(title=rf'$\alpha={a:g}$',xlim=(-5,5),ylim=(-1.1,1.1))
    if a==1:assert np.allclose(y,-x*f,atol=1e-11)
    if a==2:assert np.allclose(y,(x*x-1)*f,atol=1e-10)
save(fig,'chapter-discrete-fourier-transform','fractional-differentiation')
