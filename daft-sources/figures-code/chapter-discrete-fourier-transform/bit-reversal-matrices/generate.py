from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "listings"))
from figure_style import *
from fft_algorithms import bit_reversal_indices
fig,axes=plt.subplots(1,4,figsize=(10,2.7))
for ax,n in zip(axes,[32,64,128,256]):
    p=bit_reversal_indices(n)
    assert np.array_equal(p[p],np.arange(n))
    ax.scatter(p,np.arange(n),s=3 if n>64 else 6,color=TEAL,linewidths=0)
    ax.set(xlim=(-.5,n-.5),ylim=(n-.5,-.5),aspect='equal',title=f'$N={n}$',xticks=[0,n-1],yticks=[0,n-1])
    axes_box(ax,'column $j$','row $i$')
save(fig,'chapter-discrete-fourier-transform','bit-reversal-matrices')
