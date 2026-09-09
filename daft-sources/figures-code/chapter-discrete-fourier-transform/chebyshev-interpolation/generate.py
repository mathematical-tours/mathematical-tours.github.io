from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from scipy.interpolate import BarycentricInterpolator
from spectral_methods import chebyshev_interpolate
n=11;grid=np.linspace(-1,1,1001);f=lambda x:1/(.3**2+x*x)
fig,axes=plt.subplots(1,2,figsize=(8.7,3.1),sharey=True)
for ax,nodes,title in zip(axes,[np.linspace(-1,1,n),np.cos((np.arange(n)+.5)*np.pi/n)],['Equally spaced nodes','Roots of $T_{11}$']):
    values=f(nodes);interpolant=BarycentricInterpolator(nodes,values,rng=0)(grid)
    if title.startswith('Roots'):assert np.allclose(interpolant,chebyshev_interpolate(values,grid))
    ax.plot(grid,f(grid),'--',color=MUTED,label='$f$');ax.plot(grid,interpolant,color=TEAL,label='interpolant');ax.scatter(nodes,values,color=RED,s=16,zorder=3)
    axes_box(ax,'$x$');ax.set(title=title,xlim=(-1,1),ylim=(-.7,12.5))
axes[0].legend(frameon=False,fontsize=9)
save(fig,'chapter-discrete-fourier-transform','chebyshev-interpolation')
