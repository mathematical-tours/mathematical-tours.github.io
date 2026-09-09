from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import flower
from pde import heat_theta_2d
n=128;dt=.0002;s=dt*n*n;f=flower(n)
fig,axes=plt.subplots(3,3,figsize=(8.1,7.8))
for row,theta in enumerate([0,.5,1]):
    for col,steps in enumerate([1,2,3]):
        y=heat_theta_2d(f,s,theta,steps);ax=axes[row,col]
        ax.imshow(y,cmap='gray',vmin=0,vmax=1);axes_box(ax);ax.set(xticks=[],yticks=[])
        ax.set_xlabel(f'range [{y.min():.2g}, {y.max():.2g}]',fontsize=9)
        if row==0:ax.set_title(rf'$t={steps*dt:.4f}$')
        if col==0:ax.set_ylabel([r'Explicit: $\theta=0$',r'Crank–Nicolson: $\theta=1/2$',r'Backward Euler: $\theta=1$'][row])
        if theta>=.5:assert np.linalg.norm(y)<=np.linalg.norm(f)+1e-12
save(fig,'chapter-fourier-applications','heat-time-stepping')
