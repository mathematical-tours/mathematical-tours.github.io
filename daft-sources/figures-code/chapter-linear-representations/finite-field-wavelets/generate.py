from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools'))
sys.path.insert(0,str(ROOT/'listings'))
from figure_style import *
from representations import finite_wavelet,inverse_finite_wavelet
p=53;x=np.arange(p);center=np.arange(-26,27)
f1=np.zeros(p);f1[24]=1
f2=np.r_[np.ones(8),np.zeros(p-16),-np.ones(8)]
f3=np.cos(6*np.pi*x/p)
p1=np.zeros(p);p1[[0,1,-1]]=[2,-1,-1]
p2=np.zeros(p);p2[[0,1,-1,-2]]=[1,1,-1,-1]
p3=np.zeros(p);p3[[0,1,-1,2,-2]]=[2,-1.2,-1.2,.2,.2]
fig,axes=plt.subplots(3,3,figsize=(8.5,7.1))
for row,(f,psi) in enumerate(zip([f1,f2,f3],[p1,p2,p3])):
    W=finite_wavelet(f,psi)
    np.testing.assert_allclose(inverse_finite_wavelet(W,psi),f-f.mean(),atol=1e-13)
    a,b,c=axes[row]
    a.plot(x,f,'o-',ms=2.4,lw=.7,color=TEAL)
    discrete(b,center,np.fft.fftshift(psi),color=RED)
    a.set(xlim=(-1,53),xticks=[0,13,26,39,52]);b.set(xlim=(-27,27),xticks=[-26,-13,0,13,26])
    for ax in [a,b]:
        ax.set_box_aspect(1);axes_box(ax,'residue $x$' if ax is a else 'centered residue $x$')
        ax.margins(y=.1)
    im=c.imshow(abs(W),cmap='gray_r',origin='upper',interpolation='nearest',vmin=0,
                extent=(-.5,p-.5,p-.5,.5),aspect='equal')
    c.set(xticks=[0,13,26,39,52],yticks=[1,13,26,39,52]);axes_box(c,'translation $b$','dilation $a$')
    colorbar=fig.colorbar(im,ax=c,shrink=.8,pad=.04,fraction=.045)
    colorbar.ax.tick_params(labelsize=8,length=2)
    if row==0:
        for ax,title in zip([a,b,c],['Signal $f$',r'Wavelet $\psi$',r'$|\mathcal{W}_\psi f(b,a)|$']):ax.set_title(title)
save(fig,'chapter-linear-representations','finite-field-wavelets')
