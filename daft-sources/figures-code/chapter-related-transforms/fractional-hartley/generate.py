from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from hartley import fractional_hartley
n=32;f=np.zeros(n);f[:3]=1;k=np.fft.fftshift(np.fft.fftfreq(n)*n)
fig,axes=plt.subplots(1,5,figsize=(11,2.65),sharey=True)
for ax,t in zip(axes,np.linspace(0,1,5)):
    y=fractional_hartley(f,t);np.testing.assert_allclose(np.linalg.norm(y),np.linalg.norm(f),atol=1e-14)
    ax.plot(k,np.fft.fftshift(y.real),'o-',color=TEAL,ms=2,lw=1,label='real')
    ax.plot(k,np.fft.fftshift(y.imag),color=RED,ls='--',lw=1,label='imaginary')
    ax.set(title=rf'$t={t:g}$',xlim=(-16,15),ylim=(-.38,1.08));axes_box(ax,'index')
axes[0].set_ylabel('amplitude');axes[0].legend(fontsize=8,loc='upper left')
save(fig,'chapter-related-transforms','fractional-hartley')
