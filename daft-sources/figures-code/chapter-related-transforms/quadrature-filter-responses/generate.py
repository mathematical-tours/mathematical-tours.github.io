from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
omega=np.linspace(.012,2*np.pi-.012,2400);q=np.exp(-1j*omega)
H=[q/(1-q),(1+q)/(2*(1-q)),(1+4*q+q*q)/(3*(1-q*q))]
fig,axes=plt.subplots(1,3,figsize=(9,3.2),sharey=True)
for i,(ax,h,title) in enumerate(zip(axes,H,['Left rectangle','Trapezoidal','Simpson'])):
    if i==2:
        # Break the curve across the excluded pole, rather than joining its two sides.
        h[abs(omega-np.pi)<.009]=np.nan
    ax.plot(omega/np.pi,abs(h),color=TEAL);ax.set(title=title,xlim=(0,2),ylim=(0,3.3),xticks=[0,.5,1,1.5,2])
    axes_box(ax,r'$\omega/\pi$');ax.axvline(0,color=RED,ls=':',lw=.8);ax.axvline(2,color=RED,ls=':',lw=.8)
    if i==2:ax.axvline(1,color=RED,ls=':',lw=.8)
axes[0].set_ylabel('gain (upper range clipped)')
save(fig,'chapter-related-transforms','quadrature-filter-responses')
