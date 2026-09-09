from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from transforms import haar
n=32;t=np.linspace(0,1,n)
fig,axs=plt.subplots(2,2,figsize=(7.1,3.65))
for row,x in enumerate([np.arange(n),np.sin(2*np.pi*t)]):
    y=haar(x).real
    assert np.isclose(np.sum(x*x),np.sum(y*y))
    axs[row,0].plot(t,x,'.-',color=TEAL,markersize=4)
    discrete(axs[row,1],np.arange(n),y,BLUE)
    ymax=max(np.max(abs(y)),1)
    for boundary,color in zip([16,24,28,30,31],[TEAL,BLUE,RED,'#ac8734','#835d9e']):
        axs[row,1].axvline(boundary-.5,color=color,lw=.8,alpha=.8)
    axs[row,1].set(xlim=(-.5,31.5),ylim=(min(y)-.12*ymax,max(y)+.12*ymax),xticks=[0,8,16,24,31])
    axes_box(axs[row,0],'$t$');axes_box(axs[row,1],'Packed coefficient index')
axs[0,0].set_title('Input signal');axs[0,1].set_title('Haar coefficients')
save(fig,'chapter-duality','haar-examples')
