from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from scipy.interpolate import CubicSpline
x=np.arange(8);y=np.array([.1,1,.5,1.3,.2,.85,.8,1.3]);grid=np.linspace(0,7,701)
fig,axes=plt.subplots(1,2,figsize=(9,3.1),sharey=True)
for ax,condition,title in zip(axes,['natural','not-a-knot'],['Natural boundary conditions','Not-a-knot boundary conditions']):
    spline=CubicSpline(x,y,bc_type=condition)
    ax.plot(grid,spline(grid),color=TEAL);ax.scatter(x,y,color=RED,s=20,zorder=3)
    axes_box(ax,'$x$');ax.set(title=title,xlim=(-.15,7.15),ylim=(-.1,1.6))
    assert np.allclose(spline(x),y)
    if condition=='natural':assert np.allclose(spline([0,7],2),0)
save(fig,'chapter-fourier-applications','cubic-spline-boundaries')
