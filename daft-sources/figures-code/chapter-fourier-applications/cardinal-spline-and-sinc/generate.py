from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from splines import centered_bspline,cardinal_cubic
x=np.linspace(-5,5,2001)
fig,ax=plt.subplots(figsize=(8.2,3.3))
ax.plot(x,cardinal_cubic(x),color=TEAL,label=r'cubic cardinal kernel $\beta^3_{\rm card}$')
ax.plot(x,centered_bspline(3,x),':',color=BLUE,label=r'cubic B-spline $\beta^3$')
ax.plot(x,np.sinc(x),'--',color=RED,label=r'$\mathrm{sinc}$')
ax.plot(np.arange(-5,6),(np.arange(-5,6)==0).astype(float),'o',color=INK,ms=3)
axes_box(ax,'$x$');ax.axhline(0,color=MUTED,lw=.5);ax.set(xlim=(-5,5),ylim=(-.26,1.06));ax.legend(frameon=False,fontsize=9,loc='upper right')
save(fig,'chapter-fourier-applications','cardinal-spline-and-sinc')
