from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from chirp import scaled_fourier
n=128;x=np.zeros(n);x[:16]=1
fig,axes=plt.subplots(3,3,figsize=(9,6.5));axes=axes.ravel()
axes[0].step(np.arange(n),x,where='mid',color=TEAL);axes[0].set(title='Input signal',ylim=(-.06,1.06))
for ax,alpha in zip(axes[1:],np.arange(.6,2.01,.2)):
    y=scaled_fourier(x,alpha)
    ax.plot(y.real,color=TEAL,lw=1.1,label='real');ax.plot(y.imag,color=RED,ls='--',lw=.8,label='imaginary')
    ax.set(title=rf'$\alpha={alpha:.1f}$',ylim=(-12,17));axes_box(ax,'output index $k$')
for ax in axes:axes_box(ax);ax.set_xlim(0,n-1)
axes[0].set_xlabel('sample index $n$');axes[1].legend(fontsize=8,loc='upper right')
save(fig,'chapter-related-transforms','scaled-fourier-signal')
