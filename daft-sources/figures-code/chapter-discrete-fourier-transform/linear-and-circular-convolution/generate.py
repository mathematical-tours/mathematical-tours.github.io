from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
n=16;t=np.arange(n);f=np.minimum(t,n-t).astype(float)
k=np.arange(-2,3);g=np.ones(5)/5
linear=np.convolve(f,g)[2:2+n]
circular=sum(gj*np.roll(f,int(kj)) for kj,gj in zip(k,g))
assert np.allclose(linear[2:-2],circular[2:-2])
fig,axes=plt.subplots(2,2,figsize=(8.5,4.5))
for ax,x,y,title in zip(axes.flat,[t,k,t,t],[f,g,circular,linear],['Input $f$','Box filter $g$','Circular convolution','Zero-extended linear convolution']):
    discrete(ax,x,y);ax.set_title(title);axes_box(ax,'sample index');ax.set_ylim(bottom=0)
for ax in [axes[0,0],axes[1,0],axes[1,1]]: ax.set(ylim=(0,8.6),xticks=[0,4,8,12,15])
save(fig,'chapter-discrete-fourier-transform','linear-and-circular-convolution')
