from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"tools"))
sys.path.insert(0,str(ROOT/"listings"))
from figure_style import *
t=np.linspace(0,2*np.pi,501);zeros=np.exp(1j*np.array([-1,1])*np.pi/4);poles=.9*zeros
fig,ax=plt.subplots(figsize=(4.8,4.1));ax.plot(np.cos(t),np.sin(t),color=MUTED,lw=1)
ax.axhline(0,color=MUTED,lw=.45);ax.axvline(0,color=MUTED,lw=.45)
ax.plot(zeros.real,zeros.imag,'o',mfc='white',mec=TEAL,ms=7,mew=1.6,label='zeros')
ax.plot(poles.real,poles.imag,'x',color=RED,ms=7,mew=1.6,label='poles')
for sign in [-1,1]:
    ax.plot([0,zeros[0 if sign<0 else 1].real],[0,zeros[0 if sign<0 else 1].imag],color=MUTED,ls=':',lw=.6)
ax.set(aspect='equal',xlim=(-1.18,1.18),ylim=(-1.18,1.18));axes_box(ax,r'$\operatorname{Re}z$',r'$\operatorname{Im}z$');ax.legend(loc='upper left',fontsize=10)
save(fig,'chapter-related-transforms','notch-poles-and-zeros')
