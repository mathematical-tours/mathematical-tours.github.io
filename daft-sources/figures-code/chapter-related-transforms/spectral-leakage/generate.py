from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
n=20;nu=np.linspace(0,n,3001);k=np.arange(n)
fig,axes=plt.subplots(1,2,figsize=(9,3.3))
for ax,beta in zip(axes,[4,4.6]):
    f=np.exp(2j*np.pi*beta*k/n)
    curve=np.exp(-2j*np.pi*np.outer(nu,k)/n)@f
    ax.plot(nu,abs(curve),color=TEAL,label='finite-record spectrum')
    ax.plot(k,abs(np.fft.fft(f)),'o',ms=4,color=RED,label='DFT samples')
    ax.axvline(beta,color=MUTED,ls=':',lw=.9)
    ax.set(title=rf'$\beta={beta:g}$ cycles',xlim=(0,n),ylim=(-.5,21));axes_box(ax,r'frequency $\nu$','magnitude')
axes[0].legend(fontsize=8)
save(fig,'chapter-related-transforms','spectral-leakage')
