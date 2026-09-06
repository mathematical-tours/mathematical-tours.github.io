from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


x=np.linspace(-5,5,1800);fig,ax=plt.subplots(figsize=(10,3.2),layout='constrained');ax.set_axis_off()
for j in range(19):
    y=.032*j+.22*np.exp(-.12*(x-.07*j)**2)*np.sin(2.8*x-.26*j)
    ax.plot(x,y,color=plt.get_cmap('viridis')(.2+.5*j/18),lw=1.4,alpha=.9)
ax.set(xlim=(-5,5),ylim=(-.32,.91))
finish(fig,__file__,'A reproducible wave motif for the book cover: phase-shifted Gaussian-windowed sinusoids, with a restrained teal-to-green progression. This is a proposed cover graphic; the reading edition retains its current cover until approval.',parameters={'curves':19,'carrier_angular_frequency':2.8,'gaussian_exponent':.12})
