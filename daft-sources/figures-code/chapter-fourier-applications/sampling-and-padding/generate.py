from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
fig,axes=plt.subplots(3,3,figsize=(9,6.5),sharey=True)
for row,n in enumerate([8,16,32]):
    x=np.arange(n)/n;f=((x>=.5)&(x<1)).astype(float)
    for col,padding in enumerate([1,2,10]):
        m=n*padding;omega=2*np.pi*np.fft.fftshift(np.fft.fftfreq(m,d=1/n));F=abs(np.fft.fftshift(np.fft.fft(f,m)))/n
        ax=axes[row,col];grid=np.linspace(-np.pi*n,np.pi*n,2001)
        ax.plot(grid,.5*abs(np.sinc(grid/(4*np.pi))),color=MUTED,lw=.8,alpha=.65)
        ax.plot(omega,F,'o-',color=TEAL,ms=2,lw=.8)
        axes_box(ax,r'angular frequency $\xi$');ax.set(title=rf'$N={n},\ M={m}$',xlim=(-np.pi*n,np.pi*n),ylim=(-.025,.55))
        if col==0:ax.set_ylabel('magnitude')
save(fig,'chapter-fourier-applications','sampling-and-padding')
