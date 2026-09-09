from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from related_transform_data import triangle
f=triangle();n=f.size;m=1024;F=np.fft.fftshift(np.fft.fft(f,m));nu=np.fft.fftshift(np.fft.fftfreq(m)*n)
fig,axes=plt.subplots(3,3,figsize=(9,6),sharex=True,sharey=True)
labels=['0',r'\pi/4',r'\pi/2',r'3\pi/4',r'\pi',r'5\pi/4',r'3\pi/2',r'7\pi/4',r'2\pi']
for ax,lam,label in zip(axes.ravel(),np.linspace(0,2*np.pi,9),labels):
    y=np.cos(lam)*F.real+np.sin(lam)*F.imag
    ax.plot(nu,y,color=TEAL);ax.axhline(0,color=MUTED,lw=.5)
    ax.set(title=rf'$\lambda={label}$',xlim=(-8,8),ylim=(-70,70));axes_box(ax,r'frequency $\nu$')
for ax in axes[:,0]:ax.set_ylabel('amplitude')
save(fig,'chapter-related-transforms','generalized-hartley')
