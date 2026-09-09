from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
k=np.arange(-25,26);x=k/25
kernels=[np.exp(-10*x*x),np.cos(np.pi*x/2)**2,np.cos(np.pi*x/2)]
labels=[r'Gaussian',r'Cosine squared',r'Cosine']
omega=np.linspace(-.3*np.pi,.3*np.pi,1501)
fig,axes=plt.subplots(3,3,figsize=(9,6.4))
for row,(g,label) in enumerate(zip(kernels,labels)):
    g/=g.sum();G=np.cos(omega[:,None]*k)@g
    discrete(axes[row,0],k,g);axes[row,0].set_ylabel(label)
    axes[row,1].plot(omega/np.pi,G,color=TEAL);axes[row,1].axhline(0,color=MUTED,lw=.5)
    axes[row,2].plot(k,np.cumsum(g),color=TEAL)
    for j,ax in enumerate(axes[row]):axes_box(ax,'$n$' if j!=1 else r'$\omega/\pi$')
    axes[row,1].set_ylim(-.12,1.06);axes[row,2].set_ylim(-.03,1.03)
for ax,title in zip(axes[0],['Impulse response','Frequency response','Step response']):ax.set_title(title)
save(fig,'chapter-fourier-applications','filter-responses')
