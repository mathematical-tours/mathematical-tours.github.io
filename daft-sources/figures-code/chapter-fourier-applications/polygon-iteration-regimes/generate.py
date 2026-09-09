from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import polygon,polygon_evolution,plot_polygon
z=polygon();n=len(z);steps=np.arange(26)
F=[np.full(n,.8,dtype=complex) for _ in range(4)]
F[0][:4]=[1.05*np.exp(1j*np.pi/15),1.1*np.exp(1j*np.pi/10),1.05*np.exp(1j*np.pi/8),1.15*np.exp(1j*np.pi/12)]
F[1][:4]=[1.2*np.exp(1j*np.pi/10),1.1*np.exp(1j*np.pi/8),1.15*np.exp(1j*np.pi/4),1.15*np.exp(1j*np.pi/6)]
F[2][:3]=[1,.9*np.exp(1j*np.pi/20),.9]
F[3][:3]=[np.exp(1j*np.pi/15),np.exp(1j*np.pi/8),np.exp(1j*np.pi/6)]
fig,axes=plt.subplots(2,2,figsize=(8.5,7))
for ax,multiplier,title in zip(axes.flat,F,['(a) Expanding modes','(b) Expansion and drift','(c) Convergence','(d) Bounded rotation']):
    trajectory=polygon_evolution(z,multiplier,steps)
    for i,poly in enumerate(trajectory):plot_polygon(ax,poly,color=plt.cm.viridis(.1+.75*i/25),lw=.65,alpha=.8)
    plot_polygon(ax,trajectory[-1],color=TEAL,lw=1.4)
    axes_box(ax,'real coordinate','imaginary coordinate');ax.set(title=title,aspect='equal')
axes[1,0].plot(z.mean().real,z.mean().imag,'o',color=RED,ms=4)
save(fig,'chapter-fourier-applications','polygon-iteration-regimes')
