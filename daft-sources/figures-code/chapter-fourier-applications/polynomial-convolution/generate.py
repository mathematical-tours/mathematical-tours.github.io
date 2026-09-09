from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
p=np.array([1,2,0,1,-1,1]);q=np.array([0,1,-1,2,0,2]);product=np.convolve(p,q)
assert np.allclose(np.fft.ifft(np.fft.fft(p,11)*np.fft.fft(q,11)).real,product)
fig,axes=plt.subplots(1,3,figsize=(9,3))
for ax,coefficients,title in zip(axes,[np.pad(p,(0,5)),np.pad(q,(0,5)),product],['$P$','$Q$','$PQ$']):
    discrete(ax,np.arange(11),coefficients);axes_box(ax,'degree','coefficient');ax.set(title=title,xlim=(-.5,10.5),xticks=[0,2,4,6,8,10]);ax.axhline(0,color=MUTED,lw=.5)
save(fig,'chapter-fourier-applications','polynomial-convolution')
