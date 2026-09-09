from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from related_transform_data import notch,notch_impulse
omega=np.linspace(0,np.pi,1601);h=notch_impulse(36)
fig,axes=plt.subplots(1,2,figsize=(9,3))
axes[0].plot(omega/np.pi,abs(notch(np.exp(1j*omega))),color=TEAL);axes[0].axvline(.25,color=RED,ls=':',lw=.9)
axes[0].plot([.25],[0],'o',color=RED,ms=4);axes[0].set(title='Frequency response',xlim=(0,1),ylim=(-.03,1.17))
discrete(axes[1],np.arange(h.size),h);axes[1].set(title='Impulse response',xlim=(-1,35.5),ylim=(-.22,1.08))
axes_box(axes[0],r'$\omega/\pi$',r'$|H(e^{i\omega})|$');axes_box(axes[1],'sample index $n$','$h[n]$')
assert abs(notch(np.exp(1j*np.pi/4)))<1e-13
save(fig,'chapter-related-transforms','notch-responses')
