"""Self-convolution with cyclic indexing, without an undocumented shift."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'tools'))
from figure_style import *

x = np.arange(16)
f = ((x>=5)&(x<=10)).astype(float)
g = np.fft.ifft(np.fft.fft(f)**2).real
reference = np.array([sum(f[h]*f[(k-h)%16] for h in x) for k in x])
assert np.allclose(g,reference)
assert np.argmax(g)==15 and np.isclose(g[15],6)
fig,axs=plt.subplots(1,2,figsize=(7.0,2.55))
for ax,y,title,color in zip(axs,[f,g],[r'$f=\mathbf{1}_{\{5,\ldots,10\}}$',r'$f*f$'],[TEAL,RED]):
    discrete(ax,x,y,color)
    ax.set(title=title,xlim=(-.5,15.5),xticks=[0,5,10,15],ylim=(-.08,max(y)*1.08))
    axes_box(ax,r'$x\in\mathbb{Z}/16\mathbb{Z}$')
axs[1].text(.48,.88,'Wrap-around at $0$',transform=axs[1].transAxes,fontsize=9,ha='center')
save(fig,'chapter-finite-groups','circular-convolution')
