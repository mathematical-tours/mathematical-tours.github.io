from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "listings"))
from figure_style import *
from transforms import tensor_transform,fwt
N=32;x=np.r_[np.arange(17),np.arange(15,0,-1)]
fig,axs=plt.subplots(3,3,figsize=(7.1,5.0),sharex=True,sharey=True)
for ax,m in zip(axs.ravel(),range(9)):
    alpha=m*np.pi/16
    A=np.array([[np.cos(alpha),np.sin(alpha)],[np.sin(alpha),-np.cos(alpha)]])
    y=tensor_transform(x,A).real
    assert np.isclose(np.linalg.norm(y),np.linalg.norm(x))
    assert np.max(np.abs(y)) < 50
    if m==4:assert np.allclose(y,fwt(x)/np.sqrt(N))
    if m==8:assert np.allclose(y,x[::-1])
    discrete(ax,np.arange(N),y,TEAL if m!=4 else RED)
    title=r'$\alpha=0$' if m==0 else (r'$\alpha=\pi/2$' if m==8 else rf'$\alpha={m}\pi/16$')
    if m==4:title=r'$\alpha=\pi/4$ (Walsh)'
    ax.set(title=title,xlim=(-.5,31.5),ylim=(-50,50),xticks=[0,16,31],yticks=[-40,0,40])
    axes_box(ax)
for ax in axs[-1]:ax.set_xlabel('Index')
save(fig,'chapter-duality','intermediate-walsh-transform')
