from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from pde import poisson_source,poisson_exact,poisson_dirichlet
n=32;x=np.linspace(0,1,n+1);X,Y=np.meshgrid(x,x,indexing='ij');boundary=poisson_exact(X,Y)
u=poisson_dirichlet(poisson_source(X,Y),boundary,1/n)
residual=(u[2:,1:-1]+u[:-2,1:-1]+u[1:-1,2:]+u[1:-1,:-2]-4*u[1:-1,1:-1])*n*n-poisson_source(X[1:-1,1:-1],Y[1:-1,1:-1])
assert np.max(abs(residual))<1e-9
fig=plt.figure(figsize=(5,3.8));ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(X,Y,u,cmap='viridis',edgecolor='#20374670',linewidth=.18,antialiased=True)
ax.set(xlabel='$x$',ylabel='$y$',zlabel='$u$',xlim=(0,1),ylim=(0,1),zlim=(.45,2.85),xticks=[0,.5,1],yticks=[0,.5,1])
ax.view_init(29,-133);ax.set_box_aspect((1,1,.7));ax.tick_params(pad=1)
save(fig,'chapter-fourier-applications','poisson-reference')
