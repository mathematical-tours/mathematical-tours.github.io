from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
from related_transform_data import notch,notch_impulse
from chirp import czt
n=128;k=np.arange(n+1);theta=2*np.pi*k/n
contours=[np.exp(1j*theta),np.exp((.1+2j*np.pi)*k/n)]
fig,axes=plt.subplots(1,2,figsize=(9,3.8))
for contour,label,color,style in zip(contours,['unit circle','outward spiral'],[TEAL,RED],['-','--']):
    axes[0].plot(contour.real,contour.imag,color=color,ls=style,label=label)
    axes[1].plot(k/n,abs(notch(contour)),color=color,ls=style,label=label)
    # Finite impulse-response evaluation agrees with the rational expression up to the tail.
    values=czt(notch_impulse(n),contour[1],n+1)
    np.testing.assert_allclose(values,notch(contour),atol=2e-6,rtol=2e-6)
zeros=np.exp(1j*np.array([-1,1])*np.pi/4);poles=.9*zeros
axes[0].plot(zeros.real,zeros.imag,'o',mfc='white',mec=BLUE,ms=6,label='zeros')
axes[0].plot(poles.real,poles.imag,'x',color=INK,ms=6,label='poles')
axes[0].set(aspect='equal',xlim=(-1.18,1.18),ylim=(-1.18,1.18),title='Evaluation contours')
axes[1].set(xlim=(0,1),ylim=(-.025,1.17),title='Transfer-function magnitude');axes[0].legend(fontsize=8,loc='lower left');axes[1].legend(fontsize=9,loc='lower center')
axes_box(axes[0],r'$\operatorname{Re}z$',r'$\operatorname{Im}z$');axes_box(axes[1],r'angle / $2\pi$',r'$|H(z)|$')
save(fig,'chapter-related-transforms','z-transform-contours')
