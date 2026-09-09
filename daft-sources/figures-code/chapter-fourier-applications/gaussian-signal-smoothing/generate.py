from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import gaussian_1d,irregular_signal
from scipy.signal import convolve
f=irregular_signal();n=len(f)
fig,axes=plt.subplots(2,4,figsize=(10,4.3))
for j,sigma in enumerate([0,.8,2.2,5]):
    k,g=gaussian_1d(sigma);radius=len(k)//2
    y=f.copy() if sigma==0 else convolve(np.pad(f,radius,mode='reflect'),g,mode='valid')
    discrete(axes[0,j],k,g);axes[0,j].set(title=rf'$\sigma={sigma:g}$',xlim=(-21,21),ylim=(-.03,1.06))
    axes[1,j].plot(np.arange(n),f,color=MUTED,lw=.55,alpha=.4);axes[1,j].plot(y,color=TEAL,lw=1.2)
    axes[1,j].set(xlim=(0,n-1),ylim=(-.05,2.05))
    axes_box(axes[0,j],'kernel index');axes_box(axes[1,j],'sample index')
axes[0,0].set_ylabel('$g$');axes[1,0].set_ylabel('$f*g$')
save(fig,'chapter-fourier-applications','gaussian-signal-smoothing')
