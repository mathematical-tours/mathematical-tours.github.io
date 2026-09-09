from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from splines import centered_bspline
x=np.linspace(-2.4,2.4,1201)
fig,axes=plt.subplots(1,4,figsize=(9.5,2.7),sharey=True)
for degree,ax in enumerate(axes):
    ax.plot(x,centered_bspline(degree,x),color=TEAL);axes_box(ax,'$x$');ax.set(title=rf'$\beta^{degree}$',xlim=(-2.4,2.4),ylim=(-.04,1.08),xticks=[-2,0,2])
save(fig,'chapter-fourier-applications','centered-bsplines')
