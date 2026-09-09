from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from applications_data import polygon,polygon_evolution,plot_polygon
z=polygon();n=len(z);g=np.zeros(n);g[:2]=.5;mean_multiplier=np.fft.fft(g)
one_mode=np.full(n,.8);one_mode[1]=1
fig,axes=plt.subplots(1,2,figsize=(8.5,4))
for ax,multiplier,title in zip(axes,[mean_multiplier,one_mode],['Midpoint averaging','Keeping frequency one']):
    for i,poly in enumerate(polygon_evolution(z,multiplier,range(0,61,3))):plot_polygon(ax,poly,color=plt.cm.viridis(.1+.75*i/20),lw=.7)
    axes_box(ax,'real coordinate','imaginary coordinate');ax.set(title=title,aspect='equal')
axes[0].plot(z.mean().real,z.mean().imag,'o',color=RED,ms=5,label='centroid');axes[0].legend(frameon=False,fontsize=9)
limit=np.fft.fft(z)[1]/n*np.exp(2j*np.pi*np.arange(n)/n)
plot_polygon(axes[1],limit,color=RED,lw=1.6,label='limit');axes[1].legend(frameon=False,fontsize=9)
save(fig,'chapter-fourier-applications','polygon-filtering')
